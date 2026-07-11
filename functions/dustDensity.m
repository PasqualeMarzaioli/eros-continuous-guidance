function q = dustDensity(radius, p)
%DUSTDENSITY  Nondimensional radial dust density (two Lorentzian peaks).
%
%   q(rho) = k1/(k2+(rho-rhoA)^2) + k3/(k4+(rho-rhoB)^2); the running cost
%   in the optimal-control problem is this scalar field along the path.
%
%   Author: Pasquale Marzaioli

q = p.k1 ./ (p.k2 + (radius - p.rhoA).^2) ...
    + p.k3 ./ (p.k4 + (radius - p.rhoB).^2);
end
