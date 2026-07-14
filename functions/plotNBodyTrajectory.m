%PLOTNBODYTRAJECTORY  Ecliptic and side views of the perturbed inclined orbit.
%
%   Compares n-body and Kepler paths in the xy plane and n-body yz profile.
%
%   Author: Pasquale Marzaioli

function plotNBodyTrajectory(number, nBody, kepler)
figure(number); clf;
set(gcf, 'Position', [100, 100, 592, 480]);
layout = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout, 1, [3, 1]);
plot(kepler(:, 1), kepler(:, 2), 'Color', [0.65, 0.65, 0.65], ...
    'LineStyle', '--', 'LineWidth', 1.0, ...
    'DisplayName', 'Kepler reference'); hold on;
plot(nBody(:, 1), nBody(:, 2), 'k-', 'LineWidth', 1.2, ...
    'DisplayName', 'Third-body model');
axis equal; grid on; box on;
xlabel('x / km'); ylabel('y / km');
xlim([-60, 60]); ylim([-60, 60]);
legend('Location', 'northwest');
nexttile(layout, 4);
plot(kepler(:, 2), kepler(:, 3), 'Color', [0.65, 0.65, 0.65], ...
    'LineStyle', '--', 'LineWidth', 1.0); hold on;
plot(nBody(:, 2), nBody(:, 3), 'k-', 'LineWidth', 1.2);
grid on; box on;
xlabel('y / km'); ylabel('z / km');
xlim([-60, 60]);
ylim([-5, 5]);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
end
