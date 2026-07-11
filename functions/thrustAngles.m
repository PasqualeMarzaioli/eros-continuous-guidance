function [inPlane, outOfPlane] = thrustAngles(canonical)
%THRUSTANGLES  Resolve primer control into in-plane and out-of-plane angles.
%
%   alpha = -lambda_v / ||lambda_v|| expressed in the local RTN/NTC basis.
%
%   Author: Pasquale Marzaioli

sampleCount = size(canonical, 1);
inPlane = zeros(sampleCount, 1);
outOfPlane = zeros(sampleCount, 1);

for index = 1:sampleCount
    r = canonical(index, 1:3).';
    v = canonical(index, 4:6).';
    alpha = -canonical(index, 11:13).' ...
        / norm(canonical(index, 11:13));
    radial = r / norm(r);
    normal = cross(r, v); normal = normal / norm(normal);
    transverse = cross(normal, radial);
    inPlane(index) = atan2d(dot(alpha, radial), ...
        dot(alpha, -transverse));
    outOfPlane(index) = asind(max(-1, min(1, dot(alpha, normal))));
end
end
