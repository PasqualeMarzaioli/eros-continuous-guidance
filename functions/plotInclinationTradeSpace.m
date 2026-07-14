%PLOTINCLINATIONTRADESPACE  Trade metrics versus commanded plane-change inclination.
%
%   Exposure, TOF, propellant, peak density, and mean cross-track thrust
%   as functions of Delta-i for the continued family of solutions.
%
%   Author: Pasquale Marzaioli

function metrics = plotInclinationTradeSpace(number, inclinations, solutions, ...
        times, trajectories, p)
sampleCount = numel(inclinations);
metrics = struct('exposure', zeros(1, sampleCount), ...
    'tofMinutes', zeros(1, sampleCount), ...
    'finalMass', zeros(1, sampleCount), ...
    'peakDensity', zeros(1, sampleCount), ...
    'meanCrossThrust', zeros(1, sampleCount));

for index = 1:sampleCount
    time = times{index};
    canonical = trajectories{index};
    density = dustDensity(vecnorm(canonical(:, 1:3), 2, 2), p);
    primer = canonical(:, 11:13);
    thrustDirection = -primer ./ vecnorm(primer, 2, 2);
    metrics.exposure(index) = trapz(time, density);
    metrics.tofMinutes(index) = solutions{index}(8) * p.timeUnit / 60;
    metrics.finalMass(index) = canonical(end, 7) * p.massUnit;
    metrics.peakDensity(index) = max(density);
    metrics.meanCrossThrust(index) = trapz(time, ...
        abs(thrustDirection(:, 3))) / time(end);
end

figure(number); clf;
set(gcf, 'Position', [100, 100, 900, 700]);
layout = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
plotTradeTile(layout, inclinations, metrics.exposure, ...
    'J = \int q d\tau');
plotTradeTile(layout, inclinations, metrics.tofMinutes, 't_f / min');
plotTradeTile(layout, inclinations, p.initialMassPhysical - metrics.finalMass, ...
    'Propellant / kg');
    plotTradeTile(layout, inclinations, metrics.peakDensity, 'max q');
plotTradeTile(layout, inclinations, metrics.meanCrossThrust, ...
    'mean |\alpha_C|');
nexttile(layout); axis off;
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
end
