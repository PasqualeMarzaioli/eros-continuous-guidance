%TERMINALERRORS  Convert nondimensional terminal residuals to km and m/s.
%
%   Scales position and velocity mismatches by the problem units for reporting.
%
%   Author: Pasquale Marzaioli

function [positionError, velocityError] = terminalErrors(finalCanonical, target, p)
positionError = norm(finalCanonical(1:3) - target(1:3)) * p.distanceUnit;
velocityError = norm(finalCanonical(4:6) - target(4:6)) ...
    * p.velocityUnit * 1000;
end
