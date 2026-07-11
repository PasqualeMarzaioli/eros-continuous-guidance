function differences = plotOsculatingDifferences(number, time, nBody, ...
        kepler, mu)
%PLOTOSCULATINGDIFFERENCES  Osculating element differences for a nearly circular inclined orbit.
%
%   Uses nonsingular diagnostics [a,e,i,Omega,u]; argument of latitude u
%   avoids singular periapsis angles on near-circular orbits.
%
%   Author: Pasquale Marzaioli

nBodyElements = orbitalDiagnostics(nBody, mu);
keplerElements = orbitalDiagnostics(kepler, mu);
differences = nBodyElements - keplerElements;
differences(:, 3:5) = rad2deg(differences(:, 3:5));
days = time / 86400;

% The two-body reference must preserve its energy-derived semimajor axis.
assert(max(abs(keplerElements(:, 1) - keplerElements(1, 1))) < 1e-6);

figure(number); clf;
set(gcf, 'Position', [100, 100, 900, 700]);
layout = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
labels = {'\Delta a / km', '\Delta e', '\Delta i / deg', ...
    '\Delta \Omega / deg', '\Delta u / deg'};
for index = 1:5
    nexttile(layout);
    plot(days, differences(:, index), 'k-', 'LineWidth', 1.1);
    grid on; box on; xlabel('t / day'); ylabel(labels{index});
end
nexttile(layout); axis off;
set(findall(gcf, 'Type', 'axes'), 'FontSize', 9);
end
