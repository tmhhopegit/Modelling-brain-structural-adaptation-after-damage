function P = load_ploras(cfg)
%LOAD_PLORAS PLORAS scans: demographics, scores and binary lesions (the training data).
%   P = load_ploras(cfg) reads cfg.ploras.spreadsheet (one row per scan) and
%   each row's lesion image. Replaces the data-loading part of PLORAS_Raw.m.
%
%   P.ids         N x 1 cellstr (patients with several scans have several rows)
%   P.info        table: time_post, age, right_handed, female, english,
%                 volume, volume_left, volume_right, first_scan, lesion_ok
%   P.scores      N x S (cfg.scores.ploras), P.score_names
%   P.lesions     N x V logical; P.voxels, P.dims, P.template
%   Rows whose image is missing or unreadable have lesion_ok = false and
%   volume NaN, so they drop out of the analyses instead of counting as
%   patients without a lesion. Needs SPM for the images.
c = cfg.ploras;
T = readtable(c.spreadsheet, 'VariableNamingRule', 'preserve');
need = {c.id, c.age, c.time_post, c.sex, c.handedness, c.language, c.lesion_file};
if ~isempty(c.order), need{end + 1} = c.order; end
need = [need, cfg.scores.ploras];
missing = need(~ismember(need, T.Properties.VariableNames));
if ~isempty(missing)
    error('lesion_growth:columns', 'Columns not found in %s: %s', c.spreadsheet, strjoin(missing, ', '));
end
P.dataset = 'PLORAS';
P.ids = cellstr(string(T.(c.id)));
P.score_names = cfg.scores.names;
P.scores = numeric_columns(T, cfg.scores.ploras);

files = cellstr(string(T.(c.lesion_file)));
files(ismember(files, {'', 'NaN', '<missing>'})) = {''};
if ~isempty(c.lesion_folder)
    has = ~cellfun(@isempty, files);
    files(has) = fullfile(c.lesion_folder, files(has));
end
[P.lesions, ok, P.voxels, P.dims, P.template] = read_lesion_images(files, cfg.images.brain_mask);
[volume, left, right] = lesion_volumes(P.lesions, P.voxels, P.dims, cfg.images.midline_x);
volume(~ok) = NaN; left(~ok) = NaN; right(~ok) = NaN;

time_post = numeric_columns(T, {c.time_post});
if ~isempty(c.order)
    first_scan = numeric_columns(T, {c.order}) == 1;
else
    first_scan = earliest_per_patient(P.ids, time_post);
end
P.info = table(time_post, numeric_columns(T, {c.age}), ...
    code_demographic(T.(c.handedness), 'right_handed'), code_demographic(T.(c.sex), 'female'), ...
    code_demographic(T.(c.language), 'english'), volume, left, right, first_scan, ok, ...
    'VariableNames', {'time_post', 'age', 'right_handed', 'female', 'english', ...
    'volume', 'volume_left', 'volume_right', 'first_scan', 'lesion_ok'});
P.lesion_files = files;
end

function X = numeric_columns(T, names)
X = nan(height(T), numel(names));
for k = 1:numel(names)
    v = T.(names{k});
    if isnumeric(v) || islogical(v)
        X(:, k) = double(v);
    else
        X(:, k) = str2double(string(v));
    end
end
end

function first = earliest_per_patient(ids, t)
first = false(numel(ids), 1);
[u, ~, g] = unique(ids);
for k = 1:numel(u)
    rows = find(g == k);
    [~, j] = min(t(rows));
    first(rows(j)) = true;
end
end
