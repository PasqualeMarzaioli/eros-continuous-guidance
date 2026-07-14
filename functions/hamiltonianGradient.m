%HAMILTONIANGRADIENT  Analytical gradient of H with respect to the canonical state.
%
%   Provides dH/d[r,v,m,lambda_r,lambda_v,lambda_m] for the free-time
%   transversality row of the shooting Jacobian.
%
%   Author: Pasquale Marzaioli

function gradient = hamiltonianGradient(canonical, p)
r = canonical(1:3);
v = canonical(4:6);
m = canonical(7);
lambdaR = canonical(8:10);
lambdaV = canonical(11:13);
radius = norm(r);
primerNorm = norm(lambdaV);
gravityMatrix = eye(3) / radius^3 - 3 * (r * r.') / radius^5;
[~, dustGradient] = dustDerivatives(radius, r, p);

gradient = [
    dustGradient - gravityMatrix * lambdaV;
    lambdaR;
    p.thrust * primerNorm / m^2;
    v;
    -r / radius^3 - p.thrust * lambdaV / (m * primerNorm);
    p.massRate];
end
