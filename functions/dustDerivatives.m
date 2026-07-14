%DUSTDERIVATIVES  Synthetic radial cost with Cartesian gradient and Hessian.
%
%   Chain-rules the radial profile q(rho) into dq/dr and d^2q/dr^2 for the
%   costate equations and the analytical canonical Jacobian.
%
%   Author: Pasquale Marzaioli

function [q, gradient, hessian] = dustDerivatives(radius, r, p)
dA = radius - p.rhoA;
dB = radius - p.rhoB;
denominatorA = p.k2 + dA^2;
denominatorB = p.k4 + dB^2;
q = p.k1 / denominatorA + p.k3 / denominatorB;
qPrime = -2 * p.k1 * dA / denominatorA^2 ...
    - 2 * p.k3 * dB / denominatorB^2;
radialDirection = r / radius;
gradient = qPrime * radialDirection;

if nargout > 2
    qSecond = -2 * p.k1 / denominatorA^2 ...
        + 8 * p.k1 * dA^2 / denominatorA^3 ...
        - 2 * p.k3 / denominatorB^2 ...
        + 8 * p.k3 * dB^2 / denominatorB^3;
    hessian = qSecond * (radialDirection * radialDirection.') ...
        + (qPrime / radius) * (eye(3) ...
        - radialDirection * radialDirection.');
end
end
