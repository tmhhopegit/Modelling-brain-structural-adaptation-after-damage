function G = learn_lesion_growth(L, cfg)
%LEARN_LESION_GROWTH One Gaussian-process model per later lesion measure.
%   G = learn_lesion_growth(L, cfg), with L from longitudinal_pairs.
%   G.models{j} predicts L.Y(:, j) from L.X (earlier lesion, demographics and
%   the later time post stroke), as the original's mdl{i}.
%   G.cv (table, if cfg.growth.cv_folds > 0): how well each measure is
%   predicted for pairs the model did not see (R^2 and correlation), and how
%   that compares with simply assuming no change (the earlier value).
%   Pairs from one patient are kept in the same fold, so a patient's other
%   scans cannot leak into their own predictions. (The original did not
%   check the growth models at all.)
n_out = size(L.Y, 2);
if size(L.X, 1) < 5
    error('lesion_growth:pairs', 'Only %d longitudinal pairs; cannot learn lesion growth.', size(L.X, 1));
end
G.models = cell(n_out, 1);
for j = 1:n_out
    G.models{j} = fit_gp(L.X, L.Y(:, j), L.categorical, cfg);
end
G.y_names = L.y_names;
G.x_names = L.x_names;
G.lesion_measures = L.lesion_measures;
G.categorical = L.categorical;
G.cv = table();
if cfg.growth.cv_folds > 0
    [~, ~, patient] = unique(L.pairs.id);
    n_pat = max(patient);
    folds = min(cfg.growth.cv_folds, n_pat);
    fold_of_patient = mod(randperm(n_pat), folds) + 1;
    fold = fold_of_patient(patient)';
    pred = nan(size(L.Y));
    for f = 1:folds
        te = fold == f;
        for j = 1:n_out
            mdl = fit_gp(L.X(~te, :), L.Y(~te, j), L.categorical, cfg);
            pred(te, j) = predict(mdl, L.X(te, :));
        end
    end
    % "no change" baseline: the earlier value of the same measure
    [~, earlier_col] = ismember(L.y_names, L.x_names);
    rows = cell(n_out, 1);
    for j = 1:n_out
        y = L.Y(:, j);
        r2 = 1 - sum((y - pred(:, j)).^2) / sum((y - mean(y)).^2);
        r2_same = NaN;
        if earlier_col(j) > 0
            e = L.X(:, earlier_col(j));
            r2_same = 1 - sum((y - e).^2) / sum((y - mean(y)).^2);
        end
        rows{j} = {L.y_names{j}, r2, corr(y, pred(:, j)), r2_same};
    end
    G.cv = cell2table(vertcat(rows{:}), 'VariableNames', {'Measure', 'R2', 'r', 'R2_if_no_change'});
    G.cv_predictions = pred;
end
end
