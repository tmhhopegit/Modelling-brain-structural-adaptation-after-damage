function C = load_chile(cfg, dims)
%LOAD_CHILE Chilean patients: demographics, scores and acute lesions (the test data).
%   C = load_chile(cfg)
%   C = load_chile(cfg, dims)   require the images to have this size (e.g. P.dims)
%
%   Replaces create_chile_mat.m. Lesion images are the files matching
%   cfg.chile.lesion_pattern in cfg.chile.lesion_folder, named by patient ID;
%   they are matched to the demographics CSV by exact ID (ignoring case).
%
%   C has the same fields as load_ploras's output (ids, info, scores,
%   score_names, lesions, voxels, dims, template), plus
%   C.info.months_to_assessment  stroke -> language assessment (CAT)
%   and C.info.time_post is stroke -> scan, in months.
%   C.dropped lists patients left out (no image, no demographics, repeats).
%   Needs SPM for the images.
if nargin < 2, dims = []; end
c = cfg.chile;
files = dir(fullfile(c.lesion_folder, c.lesion_pattern));
if isempty(files)
    error('lesion_growth:images', 'No files matching %s in %s.', c.lesion_pattern, c.lesion_folder);
end
file_ids = cell(numel(files), 1);
for f = 1:numel(files)
    [~, file_ids{f}] = fileparts(files(f).name);
end
paths = fullfile({files.folder}, {files.name})';

T = readtable(c.demographics, 'VariableNamingRule', 'preserve');
need = [{c.id, c.age, c.stroke_date, c.scan_date, c.days_to_assessment, c.sex, c.handedness, c.language}, ...
    cfg.scores.chile];
missing = need(~ismember(need, T.Properties.VariableNames));
if ~isempty(missing)
    error('lesion_growth:columns', 'Columns not found in %s: %s', c.demographics, strjoin(missing, ', '));
end
ids = cellstr(string(T.(c.id)));

% one row per patient, each with an image
[~, first] = unique(lower(ids), 'stable');
repeats = setdiff(1:numel(ids), first);
[has_image, where] = ismember(lower(ids), lower(file_ids));
keep = false(numel(ids), 1);
keep(first) = true;
keep = keep & has_image;
C.dropped = struct('no_image', {ids(~has_image)}, 'repeated_rows', {ids(repeats)}, ...
    'images_without_demographics', {file_ids(~ismember(lower(file_ids), lower(ids)))});
T = T(keep, :);
C.dataset = 'Chile';
C.ids = ids(keep);
C.lesion_files = paths(where(keep));

[C.lesions, ok, C.voxels, C.dims, C.template] = read_lesion_images(C.lesion_files, cfg.images.brain_mask, dims);
[volume, left, right] = lesion_volumes(C.lesions, C.voxels, C.dims, cfg.images.midline_x);
volume(~ok) = NaN; left(~ok) = NaN; right(~ok) = NaN;

dpm = cfg.days_per_month;
months_to_scan = days(as_datetime(T.(c.scan_date)) - as_datetime(T.(c.stroke_date))) / dpm;
months_to_assessment = numeric(T.(c.days_to_assessment)) / dpm;
C.score_names = cfg.scores.names;
C.scores = nan(height(T), numel(cfg.scores.chile));
for k = 1:numel(cfg.scores.chile)
    C.scores(:, k) = numeric(T.(cfg.scores.chile{k}));
end
C.info = table(months_to_scan, numeric(T.(c.age)), ...
    code_demographic(T.(c.handedness), 'right_handed'), code_demographic(T.(c.sex), 'female'), ...
    code_demographic(T.(c.language), 'english'), volume, left, right, true(height(T), 1), ok, ...
    months_to_assessment, ...
    'VariableNames', {'time_post', 'age', 'right_handed', 'female', 'english', ...
    'volume', 'volume_left', 'volume_right', 'first_scan', 'lesion_ok', 'months_to_assessment'});
if cfg.verbose
    fprintf('Chile: %d patients with images and demographics (%d without an image, %d repeated rows).\n', ...
        numel(C.ids), numel(C.dropped.no_image), numel(repeats));
end
end

function d = as_datetime(v)
if isdatetime(v), d = v; else, d = datetime(string(v)); end
end

function x = numeric(v)
if isnumeric(v) || islogical(v), x = double(v); else, x = str2double(string(v)); end
end
