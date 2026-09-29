function preds = cv_predict_gp(X, y, categorical, cfg, folds)
%CV_PREDICT_GP Cross-validated predictions: each row predicted by a model that did not see it.
%   preds = cv_predict_gp(X, y, categorical, cfg, folds)
%   folds: number of folds (at least 2; folds >= N gives leave-one-out).
%   Uses cvpartition with the current random stream, so set rng first for
%   repeatable folds. (fitrgp's 'KFold' + kfoldPredict does the same thing.)
n = size(X, 1);
if folds >= n
    part = cvpartition(n, 'LeaveOut');
else
    part = cvpartition(n, 'KFold', folds);
end
preds = nan(n, 1);
for k = 1:part.NumTestSets
    tr = training(part, k);
    te = test(part, k);
    mdl = fit_gp(X(tr, :), y(tr), categorical, cfg);
    preds(te) = predict(mdl, X(te, :));
end
end
