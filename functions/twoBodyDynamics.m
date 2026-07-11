function derivative = twoBodyDynamics(state, mu)
%TWOBODYDYNAMICS  Unperturbed Eros-centered Kepler reference dynamics.
%
%   Baseline a = -mu r / ||r||^3 used for n-body minus Kepler comparisons.
%
%   Author: Pasquale Marzaioli

derivative = [state(4:6); -mu * state(1:3) / norm(state(1:3))^3];
end
