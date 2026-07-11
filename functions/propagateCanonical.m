function [time, trajectory] = propagateCanonical(decision, xInitial, p, ...
        odeOptions, sampleCount)
%PROPAGATECANONICAL  Reintegrate a converged PMP extremal on a uniform time grid.
%
%   Dense canonical propagation for diagnostics, cost evaluation, and plots.
%
%   Author: Pasquale Marzaioli

timeGrid = linspace(0, decision(8), sampleCount);
[time, trajectory] = ode113(@(t, y) canonicalDynamics(y, p), ...
    timeGrid, [xInitial; decision(1:7)], odeOptions);
end
