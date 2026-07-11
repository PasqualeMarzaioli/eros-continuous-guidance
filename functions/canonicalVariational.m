function derivative = canonicalVariational(augmented, p)
%CANONICALVARIATIONAL  Canonical state plus sensitivity to the seven initial costates.
%
%   Augments the 14-D PMP ODE with dS/dt = J_can * S where S is the 14-by-7
%   sensitivity of the canonical state to the free initial costates.
%
%   Author: Pasquale Marzaioli

canonical = augmented(1:14);
sensitivity = reshape(augmented(15:end), 14, 7);
derivative = [canonicalDynamics(canonical, p); ...
    reshape(canonicalJacobian(canonical, p) * sensitivity, 98, 1)];
end
