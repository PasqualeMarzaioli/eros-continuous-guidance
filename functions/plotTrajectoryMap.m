%PLOTTRAJECTORYMAP  Overlay an extremal path on the synthetic dust-cost field.
%
%   False-color q in the reference plane with start/end markers and the
%   initial/final circular orbit radii.
%
%   Author: Pasquale Marzaioli

function plotTrajectoryMap(number, ~, canonical, p, initialRadius, finalRadius)
limits = linspace(-80, 80, 400);
[xGrid, yGrid] = meshgrid(limits, limits);
radiusGrid = hypot(xGrid, yGrid);
densityGrid = dustDensity(radiusGrid / p.distanceUnit, p);
position = canonical(:, 1:3) * p.distanceUnit;
angle = linspace(0, 2 * pi, 600);

figure(number); clf;
if number == 5
    set(gcf, 'Position', [100, 100, 686, 564]);
else
    set(gcf, 'Position', [100, 100, 654, 564]);
end
surface(xGrid, yGrid, zeros(size(xGrid)), densityGrid, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
view(2); hold on;
plot(position(:, 1), position(:, 2), 'k-', 'LineWidth', 1.5);
plot(initialRadius * cos(angle), initialRadius * sin(angle), ...
    'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(finalRadius * cos(angle), finalRadius * sin(angle), ...
    'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
startHandle = plot(position(1, 1), position(1, 2), 'ko', 'MarkerSize', 9, ...
    'LineWidth', 1.2, 'DisplayName', 'Start');
endHandle = plot(position(end, 1), position(end, 2), 'kx', 'MarkerSize', 9, ...
    'LineWidth', 1.2, 'DisplayName', 'End');
axis equal; xlim([-80, 80]); ylim([-80, 80]);
grid on; box on; colormap(parula);
colorbarHandle = colorbar;
colorbarHandle.Label.String = 'q(\rho)';
xlabel('x / km'); ylabel('y / km');
legend([startHandle, endHandle], {'Start', 'End'}, 'Location', 'northwest');
set(gca, 'FontSize', 11);
end
