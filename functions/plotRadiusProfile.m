%PLOTRADIUSPROFILE  Radial history beside the outer dust-density peak profile.
%
%   Left panel: q vs rho near the outer peak; right: rho(t) with rhoB marked.
%
%   Author: Pasquale Marzaioli

function plotRadiusProfile(number, time, canonical, p)
radiusPhysical = vecnorm(canonical(:, 1:3), 2, 2) * p.distanceUnit;
rho = linspace(52, 70, 800);
density = dustDensity(rho / p.distanceUnit, p);

figure(number); clf;
if number == 6
    set(gcf, 'Position', [100, 100, 956, 428]);
else
    set(gcf, 'Position', [100, 100, 944, 432]);
end
layout = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout, 1);
plot(density, rho, 'k-', 'LineWidth', 1.4); hold on;
yline(p.rhoBPhysical, 'k--', 'LineWidth', 1.0);
xlim([0, 1.05 * max(density)]); ylim([52, 70]);
grid on; box on;
xlabel('q(\rho)'); ylabel('\rho / km');

nexttile(layout, 2, [1, 3]);
plot(time * p.timeUnit, radiusPhysical, 'k-', 'LineWidth', 1.4); hold on;
yline(p.rhoBPhysical, 'k--', 'LineWidth', 1.0);
xlim([0, time(end) * p.timeUnit]); ylim([52, 70]);
grid on; box on;
xlabel('t / s'); ylabel('\rho / km');
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
end
