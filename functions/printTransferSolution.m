%PRINTTRANSFERSOLUTION  Print high-precision transfer time, mass, costates, and errors.
%
%   Tabulates tf, final mass, initial costate vector, and terminal errors.
%
%   Author: Pasquale Marzaioli

function printTransferSolution(solution, finalMass, positionError, velocityError, p)
fprintf('tf = %.10f min, mf = %.10f kg\n', ...
    solution(8) * p.timeUnit / 60, finalMass * p.massUnit);
fprintf('lambda0 =');
fprintf(' %+.10f', solution(1:7));
fprintf('\nPosition error = %.10e km\n', positionError);
fprintf('Velocity error = %.10e m/s\n\n', velocityError);
end
