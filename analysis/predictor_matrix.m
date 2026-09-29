function [X, names, categorical] = predictor_matrix(info, pc_scores, cfg)
%PREDICTOR_MATRIX Demographics followed by lesion components, as the models use them.
%   [X, names, categorical] = predictor_matrix(info, pc_scores, cfg)
%   info: table with the columns named in cfg.model.predictors
%   pc_scores: N x K lesion component scores
%   categorical: indices (into X's columns) of cfg.model.categorical
names = cfg.model.predictors(:)';
missing = names(~ismember(names, info.Properties.VariableNames));
if ~isempty(missing)
    error('lesion_growth:predictors', 'Unknown predictors: %s', strjoin(missing, ', '));
end
X = [info{:, names}, pc_scores];
names = [names, arrayfun(@(k) sprintf('PC%d', k), 1:size(pc_scores, 2), 'UniformOutput', false)];
categorical = find(ismember(names, cfg.model.categorical));
end
