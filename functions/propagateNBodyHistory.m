%PROPAGATENBODYHISTORY  Propagate Eros-centered motion with selected SPICE third bodies.
%
%   Integrates nBodyDynamics on a prescribed output grid for sensitivity checks.
%
%   Author: Pasquale Marzaioli

function state = propagateNBodyHistory(time, initialState, etInitial, ...
        muEros, names, gm, odeOptions)
[~, state] = ode113(@(t, x) nBodyDynamics(t, x, etInitial, ...
    muEros, names, gm), time, initialState, odeOptions);
end
