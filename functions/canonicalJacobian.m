function jacobian = canonicalJacobian(canonical, p)
%CANONICALJACOBIAN  Analytical 14-by-14 Jacobian of the canonical PMP vector field.
%
%   Used to integrate shooting sensitivities; includes primer derivative,
%   gravity costate coupling, and the radial dust Hessian.
%
%   Author: Pasquale Marzaioli

r = canonical(1:3);
m = canonical(7);
lambdaV = canonical(11:13);
radius = norm(r);
primerNorm = norm(lambdaV);
alphaDirection = lambdaV / primerNorm;
gravityMatrix = eye(3) / radius^3 - 3 * (r * r.') / radius^5;
primerDerivative = eye(3) / primerNorm - ...
    (lambdaV * lambdaV.') / primerNorm^3;

scalarProduct = dot(r, lambdaV);
gravityCostateJacobian = -3 * (...
    lambdaV * r.' + r * lambdaV.' + scalarProduct * eye(3)) / radius^5 ...
    + 15 * scalarProduct * (r * r.') / radius^7;
[~, ~, dustHessian] = dustDerivatives(radius, r, p);

jacobian = zeros(14);
jacobian(1:3, 4:6) = eye(3);
jacobian(4:6, 1:3) = -gravityMatrix;
jacobian(4:6, 7) = p.thrust * alphaDirection / m^2;
jacobian(4:6, 11:13) = -p.thrust * primerDerivative / m;
jacobian(8:10, 1:3) = gravityCostateJacobian - dustHessian;
jacobian(8:10, 11:13) = gravityMatrix;
jacobian(11:13, 8:10) = -eye(3);
jacobian(14, 7) = 2 * p.thrust * primerNorm / m^3;
jacobian(14, 11:13) = -p.thrust * alphaDirection.' / m^2;
end
