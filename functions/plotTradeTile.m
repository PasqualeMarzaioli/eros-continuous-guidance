%PLOTTRADETILE  One consistently formatted inclination trade-space subplot.
%
%   Helper that draws a single metric vs Delta-i with shared styling.
%
%   Author: Pasquale Marzaioli

function plotTradeTile(layout, inclination, value, label)
nexttile(layout);
plot(inclination, value, 'ko-', 'MarkerFaceColor', 'w', 'LineWidth', 1.2);
grid on; box on; xlabel('\Delta i / deg'); ylabel(label);
end
