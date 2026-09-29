function [P, C, cfg] = make_synthetic_lesion_data(seed)
%MAKE_SYNTHETIC_LESION_DATA Made-up PLORAS and Chile datasets in which lesions really grow.
%   [P, C, cfg] = make_synthetic_lesion_data(seed)
%   Lesions are spheres in a 12 x 12 x 12 "brain" whose radius grows by up
%   to 60% over the first months (radius r0 * (1 + 0.6 (1 - exp(-t/3)))). The
%   score depends only on the lesion volume at the time of assessment:
%   score = 80 - 0.12 * volume + noise (SD 3).
%   PLORAS: 150 patients scanned and assessed at the same time; 50 of them
%   first scanned in the first 2 months and again 3-24 months later.
%   Chile: 60 patients scanned at 0.5 months, assessed 6-36 months later,
%   when their lesions are larger than on the scan. So a model trained on
%   PLORAS is over-optimistic for Chile from the acute scans, and projecting
%   the lesions forward should remove most of that. (Checked with a Python
%   copy of the analysis: acute mean error +14 to +19 over 5 seeds,
%   projected -1.9 to +1.7.)
%   cfg: config_lesion_growth with settings suited to these data.
rng(seed);
D = 12;
[gx, gy, gz] = ndgrid(1:D, 1:D, 1:D);
grid = [gx(:), gy(:), gz(:)] - 1;
growth = @(r0, t) r0 .* (1 + 0.6 * (1 - exp(-t / 3)));
sphere = @(c, r) (sqrt(sum((grid - c).^2, 2)) <= r)';

n_p = 150; n_long = 50;
ids = {}; time_post = []; order = []; lesions = false(0, D^3);
for i = 1:n_p
    c = 4 + 4 * rand(1, 3); r0 = 1.5 + 1.5 * rand;
    if i <= n_long
        t = 0.2 + 1.8 * rand; t = [t, t + 3 + 21 * rand];
    else
        t = 0.2 + 35.8 * rand;
    end
    for o = 1:numel(t)
        ids{end + 1, 1} = sprintf('P%03d', i); %#ok<AGROW>
        time_post(end + 1, 1) = t(o); %#ok<AGROW>
        order(end + 1, 1) = o; %#ok<AGROW>
        lesions(end + 1, :) = sphere(c, growth(r0, t(o))); %#ok<AGROW>
    end
end
n = numel(ids);
volume = double(sum(lesions, 2));
P.dataset = 'PLORAS';
P.ids = ids;
P.score_names = {'Synthetic score'};
P.scores = 80 - 0.12 * volume + 3 * randn(n, 1);
P.lesions = lesions;
P.voxels = (1:D^3)';
P.dims = [D D D];
P.info = table(time_post, 40 + 40 * rand(n, 1), double(rand(n, 1) < 0.9), double(rand(n, 1) < 0.5), ...
    double(rand(n, 1) < 0.8), volume, volume / 2, volume / 2, order == 1, true(n, 1), ...
    'VariableNames', {'time_post', 'age', 'right_handed', 'female', 'english', ...
    'volume', 'volume_left', 'volume_right', 'first_scan', 'lesion_ok'});

n_c = 60;
lesions = false(n_c, D^3); final_volume = zeros(n_c, 1); months = 6 + 30 * rand(n_c, 1);
for i = 1:n_c
    c = 4 + 4 * rand(1, 3); r0 = 1.5 + 1.5 * rand;
    lesions(i, :) = sphere(c, growth(r0, 0.5));
    final_volume(i) = sum(sphere(c, growth(r0, months(i))));
end
volume = double(sum(lesions, 2));
C.dataset = 'Chile';
C.ids = arrayfun(@(i) sprintf('CHPS%03d', i), (1:n_c)', 'UniformOutput', false);
C.score_names = P.score_names;
C.scores = 80 - 0.12 * final_volume + 3 * randn(n_c, 1);
C.lesions = lesions;
C.voxels = P.voxels;
C.dims = P.dims;
C.info = table(0.5 * ones(n_c, 1), 40 + 40 * rand(n_c, 1), double(rand(n_c, 1) < 0.9), ...
    double(rand(n_c, 1) < 0.5), zeros(n_c, 1), volume, volume / 2, volume / 2, true(n_c, 1), true(n_c, 1), months, ...
    'VariableNames', {'time_post', 'age', 'right_handed', 'female', 'english', ...
    'volume', 'volume_left', 'volume_right', 'first_scan', 'lesion_ok', 'months_to_assessment'});

cfg = config_lesion_growth();
cfg.verbose = false;
cfg.seed = seed;
cfg.model.predictors = {'time_post', 'volume'};   % the synthetic score depends on nothing else
cfg.model.categorical = {};
cfg.growth.age_advances = false;
end
