function S = chile_summary(C)
%CHILE_SUMMARY Sample description per score, as create_chile_mat's summary_stats.
%   S = chile_summary(C): one row per score (patients with that score): N,
%   number female, age mean and SD, and months from stroke to assessment and
%   from scan to assessment (median and SD).
n_scores = size(C.scores, 2);
rows = cell(n_scores, 9);
for k = 1:n_scores
    s = ~isnan(C.scores(:, k));
    i = C.info(s, :);
    scan_to_assessment = i.months_to_assessment - i.time_post;
    rows(k, :) = {C.score_names{k}, nnz(s), sum(i.female == 1), mean(i.age, 'omitnan'), std(i.age, 'omitnan'), ...
        median(i.months_to_assessment, 'omitnan'), std(i.months_to_assessment, 'omitnan'), ...
        median(scan_to_assessment, 'omitnan'), std(scan_to_assessment, 'omitnan')};
end
S = cell2table(rows, 'VariableNames', {'Score', 'N', 'NumFemale', 'MeanAge', 'SDAge', ...
    'MedianMonthsStrokeToAssessment', 'SDMonthsStrokeToAssessment', ...
    'MedianMonthsScanToAssessment', 'SDMonthsScanToAssessment'});
end
