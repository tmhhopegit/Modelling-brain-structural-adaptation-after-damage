function T = compare_errors(pred_a, actual_a, pred_b, actual_b, names)
%COMPARE_ERRORS Are predictions for group B more over-optimistic than for group A?
%   T = compare_errors(pred_a, actual_a, pred_b, actual_b, names)
%   names: {name of A, name of B}. Error = predicted - actual, so a positive
%   mean error means over-optimistic predictions (if higher scores are better).
%
%   T (one row): N, mean error, mean absolute error and correlation of
%   predicted with actual scores, per group; then the ANOVA of the original,
%   error ~ actual score (continuous) + group, with the p-values of both terms,
%   and the group difference in error adjusted for the actual score.
%   The actual score is in the model because errors always correlate with it
%   (regression to the mean); the group term asks whether B's errors differ
%   from A's at the same actual score.
ok_a = ~isnan(pred_a) & ~isnan(actual_a);
ok_b = ~isnan(pred_b) & ~isnan(actual_b);
ea = pred_a(ok_a) - actual_a(ok_a);
eb = pred_b(ok_b) - actual_b(ok_b);
err = [ea; eb];
actual = [actual_a(ok_a); actual_b(ok_b)];
group = [repmat(names(1), numel(ea), 1); repmat(names(2), numel(eb), 1)];
[p, tbl] = anovan(err, {actual, group}, 'continuous', 1, 'model', 'linear', ...
    'display', 'off', 'varnames', {'actual', 'group'});
% group difference at equal actual score (regression coefficient)
D = [ones(numel(err), 1), actual, double(strcmp(group, names{2}))];
beta = D \ err;
T = table(numel(ea), mean(ea), mean(abs(ea)), corr_or_nan(pred_a(ok_a), actual_a(ok_a)), ...
    numel(eb), mean(eb), mean(abs(eb)), corr_or_nan(pred_b(ok_b), actual_b(ok_b)), ...
    beta(3), tbl{3, 6}, p(2), p(1), ...
    'VariableNames', {'N_A', 'MeanError_A', 'MAE_A', 'r_A', 'N_B', 'MeanError_B', 'MAE_B', 'r_B', ...
    'AdjustedDifference', 'F_group', 'p_group', 'p_actual'});
T.Properties.Description = sprintf('A = %s, B = %s', names{1}, names{2});
end

function r = corr_or_nan(x, y)
if numel(x) < 3, r = NaN; else, r = corr(x, y); end
end
