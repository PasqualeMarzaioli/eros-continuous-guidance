%ORBITALDIAGNOSTICS  Classical elements with argument of latitude instead of omega+nu.
%
%   Returns [a, e, i, Omega, u] from energy, eccentricity vector, and the
%   ascending-node frame; unwraps Omega and u along a history.
%
%   Author: Pasquale Marzaioli

function elements = orbitalDiagnostics(state, mu)
sampleCount = size(state, 1);
elements = zeros(sampleCount, 5);
for index = 1:sampleCount
    r = state(index, 1:3).';
    v = state(index, 4:6).';
    radius = norm(r);
    h = cross(r, v);
    hDirection = h / norm(h);
    eccentricityVector = cross(v, h) / mu - r / radius;
    energy = dot(v, v) / 2 - mu / radius;
    node = [-h(2); h(1); 0];
    nodeDirection = node / norm(node);
    transverseNode = cross(hDirection, nodeDirection);
    elements(index, :) = [-mu / (2 * energy), ...
        norm(eccentricityVector), acos(hDirection(3)), ...
        atan2(nodeDirection(2), nodeDirection(1)), ...
        atan2(dot(r, transverseNode), dot(r, nodeDirection))];
end
elements(:, 4) = unwrap(elements(:, 4));
elements(:, 5) = unwrap(elements(:, 5));
end
