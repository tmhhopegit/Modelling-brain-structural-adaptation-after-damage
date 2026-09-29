%EXAMPLE_LESION_GROWTH The whole analysis, section by section.
% Edit the paths and column names below (and in config_lesion_growth.m),
% then run each section (Ctrl+Enter).
run(fullfile(fileparts(mfilename('fullpath')), 'setup_lesion_growth.m'));

cfg = config_lesion_growth();
cfg.ploras.spreadsheet  = 'C:\path\to\PLORAS.xlsx';
cfg.chile.lesion_folder = 'Z:\PLORAS_Chile_Processed\Chile_within_two_weeks';
cfg.chile.demographics  = 'Z:\Rob_Loughnan\Chile_demographics_deigo.csv';
cfg.images.brain_mask   = 'C:\path\to\brain_mask.nii';

%% Load both datasets (needs SPM); this is the slow part
P = load_ploras(cfg);
C = load_chile(cfg, P.dims);
disp(chile_summary(C))
if ~exist(cfg.results_dir, 'dir'), mkdir(cfg.results_dir); end
save(fullfile(cfg.results_dir, 'data.mat'), 'P', 'C', '-v7.3');

%% Run the analysis
R = lesion_growth_analysis(P, C, cfg);
disp(R.summary)
disp(R.growth.cv)              % is lesion change predictable at all?

%% Figures and tables
plot_lesion_growth(R);
writetable(R.summary, fullfile(cfg.results_dir, 'summary.csv'));
writetable(R.growth.cv, fullfile(cfg.results_dir, 'growth_cv.csv'));

%% Sensitivity: components from PLORAS only, and without the language predictor
% (all Chilean patients are Spanish speakers, a value few PLORAS patients have)
cfg2 = cfg;
cfg2.model.pca_fit_on = 'ploras';
cfg2.model.predictors = setdiff(cfg.model.predictors, {'english'}, 'stable');
cfg2.model.categorical = setdiff(cfg.model.categorical, {'english'}, 'stable');
R2 = lesion_growth_analysis(P, C, cfg2);
disp(R2.summary)
