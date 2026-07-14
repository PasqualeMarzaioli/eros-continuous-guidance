%HAMILTONIAN  Autonomous minimum-principle Hamiltonian for the dust cost.
%
%   H = q(rho) + lambda_r·v + lambda_v·(-r/r^3) - T||lambda_v||/m + lambda_m * mdot.
%   For free final time, a minimizing extremal satisfies H(tf) = 0.
%
%   Author: Pasquale Marzaioli

function value = hamiltonian(canonical, p)
r = canonical(1:3);
v = canonical(4:6);
m = canonical(7);
lambdaR = canonical(8:10);
lambdaV = canonical(11:13);
lambdaM = canonical(14);
radius = norm(r);
value = dustDensity(radius, p) + dot(lambdaR, v) ...
    + dot(lambdaV, -r / radius^3) ...
    - p.thrust * norm(lambdaV) / m + lambdaM * p.massRate;
end
