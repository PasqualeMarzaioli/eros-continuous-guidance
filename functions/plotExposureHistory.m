%PLOTEXPOSUREHISTORY  Instantaneous and cumulative dust exposure for two extremals.
%
%   Compares planar and inclined optima against a kinematic monotonic-radius
%   yardstick that is not dynamically feasible.
%
%   Author: Pasquale Marzaioli

function plotExposureHistory(number, timePlanar, planar, timeInclined, ...
        inclined, p, initialRadius, finalRadius)
radiusPlanar = vecnorm(planar(:, 1:3), 2, 2);
radiusInclined = vecnorm(inclined(:, 1:3), 2, 2);
densityPlanar = dustDensity(radiusPlanar, p);
densityInclined = dustDensity(radiusInclined, p);
cumulativePlanar = cumtrapz(timePlanar, densityPlanar);
cumulativeInclined = cumtrapz(timeInclined, densityInclined);

% This same-duration monotonic-radius path is explicitly kinematic: it is an
% exposure yardstick, not a dynamically feasible replacement trajectory.
referenceRadius = linspace(initialRadius / p.distanceUnit, ...
    finalRadius / p.distanceUnit, numel(timePlanar)).';
referenceDensity = dustDensity(referenceRadius, p);
referenceCumulative = cumtrapz(timePlanar, referenceDensity);

figure(number); clf;
set(gcf, 'Position', [100, 100, 850, 620]);
layout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout);
plot(timePlanar / timePlanar(end), densityPlanar, 'k-', 'LineWidth', 1.3);
hold on;
plot(timeInclined / timeInclined(end), densityInclined, 'Color', ...
    [0.85, 0.325, 0.098], 'LineWidth', 1.2);
plot(timePlanar / timePlanar(end), referenceDensity, 'k--', 'LineWidth', 1.0);
grid on; box on; ylabel('q(\rho)');
legend({'Planar extremal', 'Inclined extremal', ...
    'Monotonic-radius yardstick'}, 'Location', 'northwest');
nexttile(layout);
plot(timePlanar / timePlanar(end), cumulativePlanar, 'k-', 'LineWidth', 1.3);
hold on;
plot(timeInclined / timeInclined(end), cumulativeInclined, 'Color', ...
    [0.85, 0.325, 0.098], 'LineWidth', 1.2);
plot(timePlanar / timePlanar(end), referenceCumulative, 'k--', 'LineWidth', 1.0);
grid on; box on; xlabel('t / t_f'); ylabel('\int_0^t q(\rho) d\tau');
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
end
