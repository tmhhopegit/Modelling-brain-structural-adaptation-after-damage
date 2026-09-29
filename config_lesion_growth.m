function cfg = config_lesion_growth()
%CONFIG_LESION_GROWTH Settings for the PLORAS -> Chile lesion-growth analysis.
%   cfg = config_lesion_growth();  then change fields as needed.
%   Everything the original scripts hard-coded (column numbers, file paths,
%   which scores, which demographics are categorical) is set here.

% ---------------------------------------------------------------- PLORAS (training data)
% One row per scan. Columns are given by header name. The original read them
% by number with importdata; for reference, it used: age = data column 4,
% time post stroke = data column 5, scan order = data column 10, sex = text
% column 97, handedness = text column 98, first language = text column 99,
% binary lesion file = text column 192, scores = Behaviours(:, 32) and (:, 53),
% i.e. data columns 49 and 70. Please replace the names below with yours.
cfg.ploras.spreadsheet   = '';
cfg.ploras.id            = 'ID';
cfg.ploras.age           = 'Age';            % age at scan (the original's conversion to age at onset was commented out)
cfg.ploras.time_post     = 'TimePost';       % months from stroke to scan/assessment
cfg.ploras.order         = 'Order';          % 1 = a patient's first scan; '' = use the earliest time_post
cfg.ploras.sex           = 'Sex';            % 'F'/'Female' -> female
cfg.ploras.handedness    = 'Handedness';     % text starting 'Right' -> right-handed (see README)
cfg.ploras.language      = 'FirstLanguage';  % text containing 'english' -> English
cfg.ploras.lesion_file   = 'LesionBinary';   % full path of each binary lesion image
cfg.ploras.lesion_folder = '';               % if set, lesion_file is relative to this folder

% ---------------------------------------------------------------- Chile (test data)
% create_chile_mat.m's inputs: a folder of CHPS*.nii lesion images named by
% patient ID, and a demographics CSV.
cfg.chile.lesion_folder = '';                % e.g. Z:\PLORAS_Chile_Processed\Chile_within_two_weeks
cfg.chile.lesion_pattern = 'CHPS*.nii';      % file name without extension = patient ID
cfg.chile.demographics  = '';                % e.g. Z:\Rob_Loughnan\Chile_demographics_deigo.csv
cfg.chile.id            = 'ID';
cfg.chile.age           = 'AgeAtStroke';
cfg.chile.stroke_date   = 'LastStrokeDate';
cfg.chile.scan_date     = 'LastScanDate';
cfg.chile.days_to_assessment = 'DaysBetweenStrokeAndCAT';
cfg.chile.sex           = 'Gender';
cfg.chile.handedness    = 'Handedness';
cfg.chile.language      = 'FirstLanguages';

% ---------------------------------------------------------------- scores
% Pairs of columns holding the same task in the two datasets. The original
% paired PLORAS Behaviours(:,32) with RepW_raw and (:,53) with SPD_raw; check
% that they are the same measure on the same scale (raw vs T-scores).
cfg.scores.names  = {'Word repetition', 'Spoken picture description'};
cfg.scores.ploras = {'RepW', 'SPD'};
cfg.scores.chile  = {'RepW_raw', 'SPD_raw'};

% ---------------------------------------------------------------- images
cfg.images.brain_mask = '';                  % image whose non-zero voxels are kept ('' = all)
cfg.images.midline_x  = 46;                  % voxels x < this count as right hemisphere (91-voxel MNI: x index 1 = right)

% ---------------------------------------------------------------- model
% Predictors, in this order, followed by the lesion's principal components.
% Available: time_post, age, right_handed, female, english, volume, volume_left, volume_right.
cfg.model.predictors  = {'time_post', 'age', 'right_handed', 'female', 'english', 'volume'};
cfg.model.categorical = {'right_handed', 'female', 'english'};
cfg.model.pca_min_explained = 1;             % keep components explaining >= 1% of variance
cfg.model.pca_fit_on  = 'all';               % 'all' = PLORAS + Chile lesions (as the original), or 'ploras'
cfg.model.kernel      = 'squaredexponential';% fitrgp KernelFunction (e.g. 'ardsquaredexponential')
cfg.model.standardize = true;
cfg.model.cv_folds    = 10;                  % cross-validated PLORAS predictions

% ---------------------------------------------------------------- lesion growth
cfg.growth.min_months_between = 0;           % pairs of scans must be further apart than this
cfg.growth.cv_folds    = 10;                 % report how well growth is predicted (0 = skip)
cfg.growth.include_reverse_pairs = false;    % the original had (and commented out) later->earlier pairs
cfg.growth.age_advances = true;              % projected patients are older by the time projected (the original kept age as it was)

% ---------------------------------------------------------------- general
cfg.days_per_month = 30.44;                  % the original used 30.4 in one place and 30.5 in another
cfg.seed        = 1;
cfg.verbose     = true;
cfg.results_dir = fullfile(pwd, 'results');
end
