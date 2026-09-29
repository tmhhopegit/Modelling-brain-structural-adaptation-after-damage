classdef LesionGrowthTest < matlab.unittest.TestCase
    %LESIONGROWTHTEST  runtests('LesionGrowthTest')
    %   Synthetic data only; the loader test needs SPM and is skipped without it.

    methods (TestClassSetup)
        function paths(~)
            setup_lesion_growth();
        end
    end

    methods (Test)
        function codingIsTheSameForBothDatasets(tc)
            h = code_demographic({'Right'; 'Right now Left'; 'Left'; 'Ambi'; 'Left now Right'; ''}, 'right_handed');
            tc.verifyEqual(h(1:5), [1; 1; 0; 0; 0]);
            tc.verifyTrue(isnan(h(6)));
            tc.verifyEqual(code_demographic({'F'; 'Female'; 'M'; 'male'}, 'female'), [1; 1; 0; 0]);
            tc.verifyEqual(code_demographic({'English'; 'ENGLISH (NL)'; 'Spanish'}, 'english'), [1; 1; 0]);
        end

        function hemisphereVolumes(tc)
            dims = [91 109 91];
            voxels = sub2ind(dims, [1; 45; 46; 91], [50; 50; 50; 50], [40; 40; 40; 40]);
            lesions = logical([1 1 1 1; 1 0 0 0; 0 0 1 1]);
            [v, l, r] = lesion_volumes(lesions, voxels, dims, 46);
            tc.verifyEqual([v, l, r], [4 2 2; 1 0 1; 2 2 0]);
        end

        function pcaMatchesBuiltIn(tc)
            rng(3);
            X = rand(40, 300) < 0.2;
            X(:, 1:30) = X(:, 1:30) | (rand(40, 1) < 0.5);        % some structure
            m = lesion_pca_fit(X, 1);
            [coeff, score, ~, ~, explained] = pca(double(X));
            k = m.kept;
            tc.verifyEqual(m.explained(1:numel(explained)), explained, 'AbsTol', 1e-6);
            tc.verifyEqual(k, nnz(explained >= 1));
            tc.verifyEqual(abs(m.scores), abs(score(:, 1:k)), 'AbsTol', 1e-4);
            tc.verifyEqual(abs(m.coeff), abs(coeff(:, 1:k)), 'AbsTol', 1e-4);
            % projecting the training rows gives their scores; new rows are centred first
            tc.verifyEqual(lesion_pca_project(m, X), m.scores, 'AbsTol', 1e-4);
            Y = rand(5, 300) < 0.2;
            tc.verifyEqual(lesion_pca_project(m, Y), (double(Y) - mean(double(X))) * m.coeff, 'AbsTol', 1e-4);
        end

        function longitudinalPairsFromUnsortedRows(tc)
            [P, ~, cfg] = make_synthetic_lesion_data(1);
            pcs = zeros(numel(P.ids), 2);
            % a patient with three scans, rows in reverse time order, spread through the table
            rows = [3; 10; 20];
            P.ids(rows) = {'X1'};
            P.info.time_post(rows) = [30; 12; 2];
            P.info.volume(rows) = [300; 200; 100];
            % one unreadable scan in another longitudinal patient
            other = find(strcmp(P.ids, 'P030'));
            P.info.lesion_ok(other(2)) = false;
            P.info.volume(other(2)) = NaN;
            L = longitudinal_pairs(P, pcs, cfg);
            x = L.pairs(strcmp(L.pairs.id, 'X1'), :);
            tc.verifyEqual(height(x), 3);                            % 3 scans -> 3 pairs
            tc.verifyTrue(all(x.months_between > 0));                % always earlier -> later
            tc.verifyEqual(sort(x.months_between), [10; 18; 28]);
            tc.verifyFalse(any(strcmp(L.pairs.id, 'P030')));         % only one usable scan left
            % X ends with the later time; Y starts with the later volume
            i = find(strcmp(L.pairs.id, 'X1') & L.pairs.months_between == 28);
            tc.verifyEqual(L.X(i, end), 30);
            tc.verifyEqual(L.Y(i, 1), 300);
            tc.verifyEqual(L.X(i, strcmp(L.x_names, 'volume')), 100);
            % every other longitudinal patient is used: 50 - P030 - the three that lost a row to X1 + X1
            tc.verifyEqual(L.n_patients, 47);
            cfg.growth.min_months_between = 15;
            L2 = longitudinal_pairs(P, pcs, cfg);
            tc.verifyTrue(all(L2.pairs.months_between > 15));
        end

        function compareErrorsFindsAPlantedBias(tc)
            rng(4);
            a = 50 + 10 * randn(200, 1); b = 50 + 10 * randn(80, 1);
            pa = a + 3 * randn(200, 1);
            pb = b + 5 + 3 * randn(80, 1);                          % over-optimistic by 5
            T = compare_errors(pa, a, pb, b, {'A', 'B'});
            tc.verifyEqual(T.AdjustedDifference, 5, 'AbsTol', 1.5);
            tc.verifyLessThan(T.p_group, 1e-6);
            tc.verifyEqual(T.N_B, 80);
            T0 = compare_errors(pa, a, b + 3 * randn(80, 1), b, {'A', 'B'});
            tc.verifyLessThan(abs(T0.AdjustedDifference), 1.5);
            % missing values are skipped, not propagated
            pb(1:5) = NaN;
            tc.verifyEqual(compare_errors(pa, a, pb, b, {'A', 'B'}).N_B, 75);
        end

        function crossValidatedPredictions(tc)
            rng(5);
            x = randn(80, 2);
            y = 3 * x(:, 1) - x(:, 2) + 0.3 * randn(80, 1);
            cfg = config_lesion_growth();
            p = cv_predict_gp(x, y, [], cfg, 5);
            tc.verifyFalse(any(isnan(p)));
            tc.verifyGreaterThan(corr(p, y), 0.9);
        end

        function analysisFindsAndRemovesTheOptimism(tc)
            [P, C, cfg] = make_synthetic_lesion_data(2);
            R = lesion_growth_analysis(P, C, cfg);
            acute = R.summary(R.summary.Chile == "acute scan", :);
            projected = R.summary(R.summary.Chile == "projected to assessment", :);
            tc.verifyLessThan(abs(acute.MeanError_A), 1.5);         % PLORAS predictions unbiased
            tc.verifyGreaterThan(acute.MeanError_B, 5);             % Chile over-optimistic from acute scans
            tc.verifyLessThan(acute.p_group, 1e-3);
            tc.verifyLessThan(abs(projected.MeanError_B), acute.MeanError_B / 3);
            tc.verifyGreaterThan(acute.r_A, 0.8);
            vol = R.growth.cv(strcmp(R.growth.cv.Measure, 'volume'), :);
            tc.verifyGreaterThan(vol.R2, 0.5);                      % growth is learnable here
            tc.verifyGreaterThan(vol.R2, vol.R2_if_no_change);      % and better than assuming none
            tc.verifyEqual(height(R.scores(1).chile), 60);
            tc.verifyTrue(all(R.chile_projected.time_post == C.info.months_to_assessment));
            f = plot_lesion_growth(R);
            close(f);
        end

        function loadersReadImagesAndMatchIDs(tc)
            tc.assumeFalse(isempty(which('spm_vol')), 'SPM not on the path');
            folder = tempname; mkdir(folder); cleanup = onCleanup(@() rmdir(folder, 's'));
            dims = [10 12 10];
            V = struct('fname', '', 'dim', dims, 'dt', [spm_type('uint8') 0], 'mat', eye(4), 'pinfo', [1; 0; 0]);
            write = @(name, idx) spm_write_vol(setfield(V, 'fname', fullfile(folder, name)), ...
                reshape(double(ismember(1:prod(dims), idx)), dims)); %#ok<SFLD>
            write('a.nii', 1:10); write('b.nii', [1:5, 700:704]);
            write('CHPS01.nii', 1:3); write('CHPS02.nii', 1:6); write('CHPS99.nii', 1:2);
            T = table({'p1'; 'p1'; 'p2'}, [60; 61; 70], [3; 20; 5], [1; 2; 1], {'F'; 'F'; 'M'}, ...
                {'Left'; 'Left'; 'Right'}, {'English'; 'English'; 'Welsh'}, ...
                {fullfile(folder, 'a.nii'); fullfile(folder, 'b.nii'); fullfile(folder, 'missing.nii')}, ...
                [40; 45; 55], [30; 32; 50], 'VariableNames', ...
                {'ID', 'Age', 'TimePost', 'Order', 'Sex', 'Handedness', 'FirstLanguage', 'LesionBinary', 'RepW', 'SPD'});
            writetable(T, fullfile(folder, 'ploras.xlsx'));
            cfg = config_lesion_growth();
            cfg.verbose = false;
            cfg.ploras.spreadsheet = fullfile(folder, 'ploras.xlsx');
            cfg.images.midline_x = 6;
            tc.verifyWarning(@() load_ploras(cfg), 'lesion_growth:images');
            warning('off', 'lesion_growth:images'); restore = onCleanup(@() warning('on', 'lesion_growth:images'));
            P = load_ploras(cfg);
            tc.verifyEqual(P.info.volume(1:2), [10; 10]);
            tc.verifyTrue(isnan(P.info.volume(3)) && ~P.info.lesion_ok(3));   % missing image: excluded, not "no lesion"
            tc.verifyEqual(P.info.right_handed, [0; 0; 1]);
            tc.verifyEqual(P.info.english, [1; 1; 0]);
            tc.verifyEqual(P.info.first_scan, [true; false; true]);
            tc.verifyEqual(P.scores, [40 30; 45 32; 55 50]);

            D = table({'CHPS02'; 'chps01'; 'CHPS03'; 'CHPS01'}, [60; 70; 80; 71], ...
                datetime({'2017-01-01'; '2017-02-01'; '2017-03-01'; '2017-02-01'}), ...
                datetime({'2017-01-11'; '2017-02-05'; '2017-03-10'; '2017-02-05'}), [365; 180; 200; 180], ...
                {'M'; 'F'; 'F'; 'F'}, {'Right'; 'Left'; 'Right'; 'Left'}, {'Spanish'; 'Spanish'; 'Spanish'; 'Spanish'}, ...
                [10; 20; 30; 21], [5; 6; 7; 6], 'VariableNames', {'ID', 'AgeAtStroke', 'LastStrokeDate', ...
                'LastScanDate', 'DaysBetweenStrokeAndCAT', 'Gender', 'Handedness', 'FirstLanguages', 'RepW_raw', 'SPD_raw'});
            writetable(D, fullfile(folder, 'chile.csv'));
            cfg.chile.lesion_folder = folder;
            cfg.chile.demographics = fullfile(folder, 'chile.csv');
            C = load_chile(cfg, P.dims);
            tc.verifyEqual(C.ids, {'CHPS02'; 'chps01'});                       % exact IDs, first row of repeats
            tc.verifyEqual(C.info.volume, [6; 3]);
            tc.verifyEqual(C.info.time_post, [10; 4] / cfg.days_per_month, 'AbsTol', 1e-9);
            tc.verifyEqual(C.info.months_to_assessment, [365; 180] / cfg.days_per_month, 'AbsTol', 1e-9);
            tc.verifyEqual(C.info.right_handed, [1; 0]);
            tc.verifyEqual(C.dropped.no_image, {'CHPS03'});
            tc.verifyEqual(C.dropped.images_without_demographics, {'CHPS99'});
            tc.verifyEqual(height(chile_summary(C)), 2);
        end
    end
end
