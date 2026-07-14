%NBODYDYNAMICS  Eros two-body gravity plus differential third-body accelerations.
%
%   a = -mu_Eros r/||r||^3 + sum GM_i ( (r_i-r)/||r_i-r||^3 - r_i/||r_i||^3 )
%   with body states from SPICE in ECLIPJ2000 relative to EROS.
%
%   Author: Pasquale Marzaioli

function derivative = nBodyDynamics(time, state, etInitial, muEros, names, gm)
r = state(1:3);
acceleration = -muEros * r / norm(r)^3;
et = etInitial + time;

for index = 1:numel(names)
    bodyState = cspice_spkezr(names{index}, et, ...
        'ECLIPJ2000', 'NONE', 'EROS');
    bodyPosition = bodyState(1:3);
    relativePosition = bodyPosition - r;
    acceleration = acceleration + gm(index) * (...
        relativePosition / norm(relativePosition)^3 ...
        - bodyPosition / norm(bodyPosition)^3);
end
derivative = [state(4:6); acceleration];
end
