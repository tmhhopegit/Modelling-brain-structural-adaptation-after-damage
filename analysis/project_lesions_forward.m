function [info_proj, pcs_proj] = project_lesions_forward(G, info, pc_scores, target_months, cfg)
%PROJECT_LESIONS_FORWARD Predict what patients' lesions will look like at a later time.
%   [info_proj, pcs_proj] = project_lesions_forward(G, info, pc_scores, target_months, cfg)
%   G: from learn_lesion_growth; info/pc_scores: the patients' current scans;
%   target_months: time post stroke to project to (e.g. the assessment).
%   info_proj: info with time_post = target_months, the lesion measures
%   replaced by their projections (volumes are kept >= 0) and, if
%   cfg.growth.age_advances, age increased by the elapsed time.
%   pcs_proj: projected lesion components.
[Xnow, xnames] = predictor_matrix(info, pc_scores, cfg);
if ~isequal([xnames, {'later_time_post'}], G.x_names)
    error('lesion_growth:growth', 'The growth models were trained with different predictors.');
end
Xq = [Xnow, target_months(:)];
Yq = nan(size(Xq, 1), numel(G.models));
for j = 1:numel(G.models)
    Yq(:, j) = predict(G.models{j}, Xq);
end
info_proj = info;
for m = 1:numel(G.lesion_measures)
    info_proj.(G.lesion_measures{m}) = max(Yq(:, m), 0);
end
if cfg.growth.age_advances
    info_proj.age = info.age + (target_months(:) - info.time_post) / 12;
end
info_proj.time_post = target_months(:);
pcs_proj = Yq(:, numel(G.lesion_measures) + 1:end);
end
