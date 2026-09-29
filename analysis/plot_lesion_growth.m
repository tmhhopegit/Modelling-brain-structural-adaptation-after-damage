function fig = plot_lesion_growth(R)
%PLOT_LESION_GROWTH Predicted against actual scores: PLORAS (cross-validated), Chile from the acute scan, Chile projected.
%   fig = plot_lesion_growth(R), R from lesion_growth_analysis.
%   One row per score. Points above the dashed line are over-optimistic
%   predictions; each panel's title gives the mean error (predicted - actual).
n = numel(R.scores);
fig = figure('Name', 'Lesion growth', 'Color', 'w', 'Position', [100 100 1100 330 * n]);
for k = 1:n
    S = R.scores(k);
    panels = {S.ploras.actual, S.ploras.predicted, 'PLORAS (cross-validated)'
              S.chile.actual,  S.chile.predicted,  'Chile, acute scan'
              S.chile.actual,  S.chile.predicted_projected, 'Chile, lesion projected to assessment'};
    lims = [min([S.ploras.actual; S.chile.actual]), max([S.ploras.actual; S.chile.actual])];
    for j = 1:3
        ax = subplot(n, 3, 3 * (k - 1) + j);
        a = panels{j, 1}; p = panels{j, 2};
        ok = ~isnan(a) & ~isnan(p);
        scatter(ax, a(ok), p(ok), 18, 'filled', 'MarkerFaceAlpha', 0.6);
        hold(ax, 'on');
        plot(ax, lims, lims, 'k--');
        hold(ax, 'off');
        axis(ax, 'square'); box(ax, 'off');
        xlabel(ax, ['Actual ' S.name]); ylabel(ax, 'Predicted');
        title(ax, sprintf('%s\nmean error %+.2f, r = %.2f', panels{j, 3}, mean(p(ok) - a(ok)), corr(a(ok), p(ok))), ...
            'FontWeight', 'normal');
    end
end
end
