%PLOTNTCERRORS  N-body minus Kepler errors in the rotating NTC frame.
%
%   Resolves position/velocity differences into radial (N), tangential (T),
%   and cross-track (C) components plus inertial norms vs time.
%
%   Author: Pasquale Marzaioli

function plotNtcErrors(number, time, nBody, kepler)
[positionNtc, velocityNtc, positionDifference, velocityDifference] = ...
    ntcDifferences(nBody, kepler);

days = time / 86400;
figure(number); clf;
set(gcf, 'Position', [100, 100, 966, 698]);
layout = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
positionLabels = {'|\Delta r_N| / km', '|\Delta r_T| / km', ...
    '|\Delta r_C| / km'};
velocityLabels = {'|\Delta v_N| / km/s', '|\Delta v_T| / km/s', ...
    '|\Delta v_C| / km/s'};

for component = 1:3
    nexttile(layout, 2 * component - 1);
    semilogy(days, max(abs(positionNtc(:, component)), 1e-12), ...
        'k-', 'LineWidth', 1.0);
    grid on; box on; ylabel(positionLabels{component}); xlim([0, days(end)]);
    nexttile(layout, 2 * component);
    semilogy(days, max(abs(velocityNtc(:, component)), 1e-12), ...
        'k-', 'LineWidth', 1.0);
    grid on; box on; ylabel(velocityLabels{component}); xlim([0, days(end)]);
end

nexttile(layout, 7);
semilogy(days, max(vecnorm(positionDifference, 2, 2), 1e-12), ...
    'k-', 'LineWidth', 1.0);
grid on; box on; xlabel('t / day'); ylabel('||\Delta r_{in}|| / km');
xlim([0, days(end)]);
nexttile(layout, 8);
semilogy(days, max(vecnorm(velocityDifference, 2, 2), 1e-12), ...
    'k-', 'LineWidth', 1.0);
grid on; box on; xlabel('t / day'); ylabel('||\Delta v_{in}|| / km/s');
xlim([0, days(end)]);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 8);
end
