function state = propagateNBodyHistory(time, initialState, etInitial, ...
        muEros, names, gm, odeOptions)
%PROPAGATENBODYHISTORY  Propagate Eros-centered motion with selected SPICE third bodies.
%
%   Integrates nBodyDynamics on a prescribed output time grid for validation.
%
%   Author: Pasquale Marzaioli

[~, state] = ode113(@(t, x) nBodyDynamics(t, x, etInitial, ...
    muEros, names, gm), time, initialState, odeOptions);
end
