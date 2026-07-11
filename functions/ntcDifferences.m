function [positionNtc, velocityNtc, positionDifference, velocityDifference] = ...
%NTCDIFFERENCES  Inertial state differences resolved in the rotating NTC frame.
%
%   Applies the transport theorem: v_NTC = R (v_in - v_ref) - omega x r_NTC
%   with omega along the orbit normal equal to ||r x v|| / ||r||^2.
%
%   Author: Pasquale Marzaioli

        ntcDifferences(comparison, reference)
% Resolve inertial state differences in the rotating reference NTC frame.
sampleCount = size(reference, 1);
positionNtc = zeros(sampleCount, 3);
velocityNtc = zeros(sampleCount, 3);
positionDifference = comparison(:, 1:3) - reference(:, 1:3);
velocityDifference = comparison(:, 4:6) - reference(:, 4:6);
for index = 1:sampleCount
    rotation = ntcRotation(reference(index, 1:6).');
    localPosition = rotation * positionDifference(index, :).';
    r = reference(index, 1:3).';
    v = reference(index, 4:6).';
    frameRate = norm(cross(r, v)) / norm(r)^2;
    localVelocity = rotation * velocityDifference(index, :).' ...
        - cross([0; 0; frameRate], localPosition);
    positionNtc(index, :) = localPosition.';
    velocityNtc(index, :) = localVelocity.';
end
end
