function L = longitudinal_pairs(P, pc_scores, cfg)
%LONGITUDINAL_PAIRS Earlier/later scan pairs from PLORAS patients scanned more than once.
%   L = longitudinal_pairs(P, pc_scores, cfg)
%   For each patient (grouped by exact ID, wherever the rows are), scans with
%   a readable lesion are sorted by time post stroke and every earlier ->
%   later pair is used (as nchoosek in the original), if the later scan is
%   more than cfg.growth.min_months_between later.
%
%   L.X: early predictors (cfg.model.predictors), early components, and the
%        later scan's time post stroke    -> what the growth models see
%   L.Y: later lesion measures (volume, and volume_left/right if they are
%        predictors) and later components -> what they predict
%   L.x_names, L.y_names, L.pairs (table: id, early_row, late_row, months_between).
%
%   The original found patients by walking through consecutive rows with the
%   same ID and the 'Order' column. Its bookkeeping skipped every other
%   patient (after closing one patient's run it forgot where the next one
%   started) and never used the last patient; it also filtered on
%   ploras_dems(found_last), a row index into the wrong list.
lesion_measures = intersect({'volume', 'volume_left', 'volume_right'}, cfg.model.predictors, 'stable');
[Xall, xnames] = predictor_matrix(P.info, pc_scores, cfg);
Yall = [P.info{:, lesion_measures}, pc_scores];
ynames = [lesion_measures, arrayfun(@(k) sprintf('PC%d', k), 1:size(pc_scores, 2), 'UniformOutput', false)];
usable = P.info.lesion_ok & P.info.volume > 0 & ~isnan(P.info.time_post) & all(~isnan(Xall), 2);

[ids, ~, g] = unique(P.ids);
X = []; Y = []; rows = zeros(0, 2); ids_out = {};
for k = 1:numel(ids)
    r = find(g == k & usable);
    if numel(r) < 2, continue; end
    [~, o] = sort(P.info.time_post(r));
    r = r(o);
    pairs = nchoosek(r, 2);                     % earlier, later (sorted by time)
    if cfg.growth.include_reverse_pairs
        pairs = [pairs; fliplr(pairs)];         %#ok<AGROW>
    end
    gap = P.info.time_post(pairs(:, 2)) - P.info.time_post(pairs(:, 1));
    pairs = pairs(abs(gap) > cfg.growth.min_months_between, :);
    if isempty(pairs), continue; end
    X = [X; Xall(pairs(:, 1), :), P.info.time_post(pairs(:, 2))]; %#ok<AGROW>
    Y = [Y; Yall(pairs(:, 2), :)];                                 %#ok<AGROW>
    rows = [rows; pairs];                                           %#ok<AGROW>
    ids_out = [ids_out; repmat(ids(k), size(pairs, 1), 1)];         %#ok<AGROW>
end
L.X = X;
L.Y = Y;
L.x_names = [xnames, {'later_time_post'}];
L.y_names = ynames;
[~, cat_idx] = ismember(cfg.model.categorical, L.x_names);
L.categorical = cat_idx(cat_idx > 0);
L.lesion_measures = lesion_measures;
L.pairs = table(ids_out, rows(:, 1), rows(:, 2), ...
    P.info.time_post(rows(:, 2)) - P.info.time_post(rows(:, 1)), ...
    'VariableNames', {'id', 'early_row', 'late_row', 'months_between'});
L.n_patients = numel(unique(ids_out));
end
