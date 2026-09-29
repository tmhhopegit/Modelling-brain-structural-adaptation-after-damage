function mdl = fit_gp(X, y, categorical, cfg)
%FIT_GP Gaussian-process regression with the configured kernel.
%   mdl = fit_gp(X, y, categorical, cfg)  (fitrgp; Statistics and Machine Learning Toolbox)
args = {'KernelFunction', cfg.model.kernel, 'Standardize', cfg.model.standardize};
if ~isempty(categorical)
    args = [args, {'CategoricalPredictors', categorical}];
end
mdl = fitrgp(X, y, args{:});
end
