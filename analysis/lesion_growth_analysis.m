function R = lesion_growth_analysis(P, C, cfg)
%LESION_GROWTH_ANALYSIS Do acute scans make predictions over-optimistic, and does projecting the lesion forward fix it?
%   R = lesion_growth_analysis(P, C, cfg)
%   P: PLORAS (load_ploras), C: Chile (load_chile), cfg: config_lesion_growth.
%   Replaces Rob_Analysis_Chilean.m:
%
%   1. Train a Gaussian process on PLORAS patients' first scans
%      (demographics + lesion components) for each score. PLORAS
%      predictions are cross-validated; Chilean patients are predicted by
%      the model trained on all PLORAS patients, from their acute scans.
%      compare_errors asks whether the Chilean predictions are more
%      over-optimistic at the same actual score.
%   2. Learn how lesions change between scans from PLORAS patients scanned
%      more than once (longitudinal_pairs, learn_lesion_growth).
%   3. Project each Chilean lesion forward to the time of the assessment,
%      predict again, and compare again.
%
%   R.summary       one row per score and condition ('acute scan',
%                   'projected to assessment'): the compare_errors table
%   R.scores(k)     per score: PLORAS ids/actual/predicted (cross-validated),
%                   Chile ids/actual/predicted/predicted_projected
%   R.pca           components kept and variance explained
%   R.growth        the pairs used, and how well growth is predicted (R.growth.cv)
%   R.chile_summary sample description (chile_summary)
%   R.models{k}     the PLORAS model for score k
if ~isequal(P.dims, C.dims) || ~isequal(P.voxels, C.voxels)
    error('lesion_growth:space', 'PLORAS and Chilean lesions are not in the same space / brain mask.');
end

% ---------------------------------------------------------------- lesion components
switch cfg.model.pca_fit_on
    case 'all'
        fit_rows = [P.lesions(P.info.lesion_ok, :); C.lesions(C.info.lesion_ok, :)];
    case 'ploras'
        fit_rows = P.lesions(P.info.lesion_ok, :);
    otherwise
        error('lesion_growth:pca', 'cfg.model.pca_fit_on must be ''all'' or ''ploras''.');
end
pca_model = lesion_pca_fit(fit_rows, cfg.model.pca_min_explained);
clear fit_rows
pcP = lesion_pca_project(pca_model, P.lesions);
pcC = lesion_pca_project(pca_model, C.lesions);
R.pca = struct('kept', pca_model.kept, 'explained', pca_model.explained(1:pca_model.kept), ...
    'fit_on', cfg.model.pca_fit_on);
say(cfg, 'Lesion components: %d, explaining %.1f%% of the variance.\n', pca_model.kept, sum(R.pca.explained));

% ---------------------------------------------------------------- 1. acute predictions
[XP, xnames, cat_idx] = predictor_matrix(P.info, pcP, cfg);
XC = predictor_matrix(C.info, pcC, cfg);
cross = P.info.first_scan & P.info.lesion_ok & P.info.volume > 0 & all(~isnan(XP), 2);
chile_ok = C.info.lesion_ok & all(~isnan(XC), 2);
n_scores = size(P.scores, 2);
R.models = cell(n_scores, 1);
summary = cell(0, 1);
for k = 1:n_scores
    s = cross & ~isnan(P.scores(:, k));
    sc = chile_ok & ~isnan(C.scores(:, k));
    say(cfg, '%s: training on %d PLORAS patients; %d Chilean patients.\n', P.score_names{k}, nnz(s), nnz(sc));
    rng(cfg.seed);
    pred_cv = cv_predict_gp(XP(s, :), P.scores(s, k), cat_idx, cfg, cfg.model.cv_folds);
    R.models{k} = fit_gp(XP(s, :), P.scores(s, k), cat_idx, cfg);
    pred_c = predict(R.models{k}, XC(sc, :));
    R.scores(k).name = P.score_names{k};
    R.scores(k).ploras = table(P.ids(s), P.scores(s, k), pred_cv, 'VariableNames', {'id', 'actual', 'predicted'});
    R.scores(k).chile = table(C.ids(sc), C.scores(sc, k), pred_c, 'VariableNames', {'id', 'actual', 'predicted'});
    T = compare_errors(pred_cv, P.scores(s, k), pred_c, C.scores(sc, k), {'PLORAS', 'Chile'});
    summary{end + 1, 1} = [table(string(P.score_names{k}), "acute scan", 'VariableNames', {'Score', 'Chile'}), T]; %#ok<AGROW>
end
R.predictor_names = xnames;

% ---------------------------------------------------------------- 2. lesion growth
L = longitudinal_pairs(P, pcP, cfg);
say(cfg, 'Lesion growth: %d pairs of scans from %d patients.\n', size(L.X, 1), L.n_patients);
rng(cfg.seed);
G = learn_lesion_growth(L, cfg);
R.growth = struct('pairs', L.pairs, 'cv', G.cv, 'models', {G.models}, 'measures', {G.y_names});
if cfg.verbose && ~isempty(G.cv), disp(G.cv); end

% ---------------------------------------------------------------- 3. projected predictions
[info_proj, pc_proj] = project_lesions_forward(G, C.info, pcC, C.info.months_to_assessment, cfg);
XCp = predictor_matrix(info_proj, pc_proj, cfg);
R.chile_projected = info_proj;
for k = 1:n_scores
    s = cross & ~isnan(P.scores(:, k));
    sc = chile_ok & ~isnan(C.scores(:, k)) & all(~isnan(XCp), 2);
    pred_p = nan(numel(C.ids), 1);
    pred_p(sc) = predict(R.models{k}, XCp(sc, :));
    [~, where] = ismember(R.scores(k).chile.id, C.ids);
    R.scores(k).chile.predicted_projected = pred_p(where);
    T = compare_errors(R.scores(k).ploras.predicted, P.scores(s, k), pred_p(sc), C.scores(sc, k), {'PLORAS', 'Chile'});
    summary{end + 1, 1} = [table(string(P.score_names{k}), "projected to assessment", ...
        'VariableNames', {'Score', 'Chile'}), T]; %#ok<AGROW>
end
R.summary = vertcat(summary{:});
R.summary = sortrows(R.summary, 'Score', 'ascend');
R.chile_summary = chile_summary(C);
if cfg.verbose
    disp(R.summary(:, {'Score', 'Chile', 'MeanError_A', 'MeanError_B', 'r_A', 'r_B', 'AdjustedDifference', 'p_group'}));
end
end

function say(cfg, varargin)
if cfg.verbose, fprintf(varargin{:}); end
end
