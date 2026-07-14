%PLOTPERTURBATIONATTRIBUTION  Signed NTC errors and dominant-source comparison.
%
%   Compares full, Sun-only, and aggregate non-Sun third-body histories.
%
%   Author: Pasquale Marzaioli

function plotPerturbationAttribution(number, time, fullState, sunState, ...
        nonSunState, kepler)
[positionFull, velocityFull] = ntcDifferences(fullState, kepler);
[positionSun, ~] = ntcDifferences(sunState, kepler);
[positionNonSun, ~] = ntcDifferences(nonSunState, kepler);
days = time / 86400;

figure(number); clf;
set(gcf, 'Position', [100, 100, 780, 850]);
layout = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
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
set(findall(gcf, 'Type', 'axes'), 'FontSize', 9);
end
