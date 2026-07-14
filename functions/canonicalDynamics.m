%CANONICALDYNAMICS  Full-thrust PMP state and costate equations for the dust problem.
%
%   Thrust direction is anti-aligned with the primer lambda_v (Lawden).
%   Costates evolve under Kepler gravity linearized about r, the dust
%   cost gradient, and the mass-flow coupling from continuous thrust.
%
%   Author: Pasquale Marzaioli

function derivative = canonicalDynamics(canonical, p)
r = canonical(1:3);
v = canonical(4:6);
m = canonical(7);
lambdaR = canonical(8:10);
lambdaV = canonical(11:13);
radius = norm(r);
primerNorm = norm(lambdaV);

if primerNorm < 1e-13 || m <= 0
    error('The canonical state reached a singular primer vector or mass.');
end

alpha = -lambdaV / primerNorm;
[~, dustGradient] = dustDerivatives(radius, r, p);
gravityMatrix = eye(3) / radius^3 - 3 * (r * r.') / radius^5;

derivative = [
    v;
    -r / radius^3 + p.thrust * alpha / m;
    p.massRate;
    gravityMatrix * lambdaV - dustGradient;
    -lambdaR;
    -p.thrust * primerNorm / m^2];
end
