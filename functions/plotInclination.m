function plotInclination(number, time, canonical, p)
%PLOTINCLINATION  Osculating inclination history from angular momentum.
%
%   i = acos(h_z / ||h||) along the converged canonical trajectory.
%
%   Author: Pasquale Marzaioli

angularMomentum = cross(canonical(:, 1:3), canonical(:, 4:6), 2);
inclination = acosd(angularMomentum(:, 3) ...
    ./ vecnorm(angularMomentum, 2, 2));
figure(number); clf;
set(gcf, 'Position', [100, 100, 594, 390]);
plot(time * p.timeUnit, inclination, 'k-', 'LineWidth', 1.4);
grid on; box on;
xlabel('t / s'); ylabel('inc / deg');
xlim([0, time(end) * p.timeUnit]); ylim([0, 3.5]);
set(gca, 'FontSize', 11);
end
