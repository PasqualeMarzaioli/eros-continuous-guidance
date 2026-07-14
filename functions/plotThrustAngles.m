%PLOTTHRUSTANGLES  Thrust steering angles in the radial-transverse-cross frame.
%
%   In-plane angle from retrograde transverse toward radial, and out-of-plane
%   elevation relative to the orbit normal, both from the primer direction.
%
%   Author: Pasquale Marzaioli

function plotThrustAngles(number, time, canonical, p)
[inPlane, outOfPlane] = thrustAngles(canonical);
figure(number); clf;
if number == 7
    set(gcf, 'Position', [100, 100, 956, 570]);
else
    set(gcf, 'Position', [100, 100, 952, 566]);
end
layout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout);
plot(time * p.timeUnit, inPlane, 'k-', 'LineWidth', 1.4);
grid on; box on;
ylabel('\alpha_{T,plane} / deg');
xlim([0, time(end) * p.timeUnit]);
nexttile(layout);
plot(time * p.timeUnit, outOfPlane, 'k-', 'LineWidth', 1.4);
grid on; box on;
xlabel('t / s'); ylabel('\alpha_{T,cross} / deg');
xlim([0, time(end) * p.timeUnit]);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
end
