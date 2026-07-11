function plotPerturbationAttribution(number, time, fullState, sunState, ...
        nonSunState, kepler, bodyNames, finalMarginalNtc)
%PLOTPERTURBATIONATTRIBUTION  Signed NTC errors and leave-one-out body contributions.
%
%   Compares full, Sun-only, and non-Sun n-body errors and bars each body's
%   marginal final tangential position contribution.
%
%   Author: Pasquale Marzaioli

[positionFull, velocityFull] = ntcDifferences(fullState, kepler);
[positionSun, ~] = ntcDifferences(sunState, kepler);
[positionNonSun, ~] = ntcDifferences(nonSunState, kepler);
days = time / 86400;

figure(number); clf;
set(gcf, 'Position', [100, 100, 1000, 720]);
layout = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout);
plot(days, positionFull, 'LineWidth', 1.0);
grid on; box on; xlabel('t / day'); ylabel('\Delta r_{NTC} / km');
legend({'N', 'T', 'C'}, 'Location', 'northwest');
nexttile(layout);
plot(days, velocityFull * 1000, 'LineWidth', 1.0);
grid on; box on; xlabel('t / day'); ylabel('\Delta v_{NTC} / m/s');
legend({'N', 'T', 'C'}, 'Location', 'northwest');
nexttile(layout);
plot(days, positionFull(:, 2), 'k-', 'LineWidth', 1.3); hold on;
plot(days, positionSun(:, 2), '--', 'Color', [0.85, 0.325, 0.098], ...
    'LineWidth', 1.2);
plot(days, positionNonSun(:, 2), 'k:', 'LineWidth', 1.2);
grid on; box on; xlabel('t / day'); ylabel('\Delta r_T / km');
legend({'All bodies', 'Sun only', 'All except Sun'}, 'Location', 'northwest');
nexttile(layout);
bar(finalMarginalNtc(:, 2), 'FaceColor', [0.35, 0.35, 0.35]);
grid on; box on; ylabel('Final marginal \Delta r_T / km');
labels = strrep(bodyNames, ' BARYCENTER', '');
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels);
xtickangle(45);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 9);
end
