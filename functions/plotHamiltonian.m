%PLOTHAMILTONIAN  Verify that the autonomous free-time Hamiltonian stays near zero.
%
%   Plots H(t) along the extremal; conservation near zero checks both the
%   integrator and satisfaction of the transversality condition.
%
%   Author: Pasquale Marzaioli

function values = plotHamiltonian(number, time, canonical, p)
values = zeros(size(time));
for index = 1:numel(time)
    values(index) = hamiltonian(canonical(index, :).', p);
end
figure(number); clf;
if number == 8
    set(gcf, 'Position', [100, 100, 600, 436]);
else
    set(gcf, 'Position', [100, 100, 552, 446]);
end
plot(time * p.timeUnit, values, 'k-', 'LineWidth', 1.4);
grid on; box on;
xlabel('t / s'); ylabel('H');
xlim([0, time(end) * p.timeUnit]);
set(gca, 'FontSize', 11);
end
