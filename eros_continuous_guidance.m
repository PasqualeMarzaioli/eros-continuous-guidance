% eros_continuous_guidance.m
% Eros continuous low-thrust transfer: PMP single shooting,
% inclination continuation, and SPICE n-body validation.
% Author: Pasquale Marzaioli

clear; close all; clc;
projectDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDirectory, 'functions'));

%% Shared physical data and nondimensionalization
p.hInitial = 52.25;          % km
p.hFinal = 36.15;            % km
p.asteroidRadius = 17.00;    % km
p.muPhysical = 3.5e-4;       % km^3/s^2
p.rhoAPhysical = 23.314 + p.asteroidRadius;
p.rhoBPhysical = 42.170 + p.asteroidRadius;
p.k1 = 7.393750e-3;
p.k2 = 7.500000e-3;
p.k3 = 3.696875e-4;
p.k4 = 6.250000e-4;
p.initialMassPhysical = 25.0; % kg
p.thrustPhysical = 21.59e-9; % kg km/s^2
p.ispPhysical = 382.82;       % s
p.g0Physical = 9.80665e-3;    % km/s^2
p.distanceUnit = p.hInitial + p.asteroidRadius;
p.massUnit = p.initialMassPhysical;
p.timeUnit = sqrt(p.distanceUnit^3 / p.muPhysical);
p.velocityUnit = p.distanceUnit / p.timeUnit;
p.mu = 1;
p.rhoA = p.rhoAPhysical / p.distanceUnit;
p.rhoB = p.rhoBPhysical / p.distanceUnit;
p.thrust = p.thrustPhysical / ...
    (p.massUnit * p.distanceUnit / p.timeUnit^2);
p.isp = p.ispPhysical / p.timeUnit;
p.g0 = p.g0Physical * p.timeUnit^2 / p.distanceUnit;
p.massRate = -p.thrust / (p.isp * p.g0);

odeOptions = odeset('RelTol', 1e-11, 'AbsTol', 1e-12);
screeningOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);

%% Dust density and circular boundary states
initialRadius = p.hInitial + p.asteroidRadius;
finalRadius = p.hFinal + p.asteroidRadius;
phase = deg2rad(45);
initialSpeed = sqrt(p.muPhysical / initialRadius);
finalSpeed = sqrt(p.muPhysical / finalRadius);

rInitialPhysical = initialRadius * [cos(phase); sin(phase); 0];
vInitialPhysical = initialSpeed * [-sin(phase); cos(phase); 0];
rFinalPhysical = [finalRadius; 0; 0];
vFinalPhysical = [0; finalSpeed; 0];

xInitial = [rInitialPhysical / p.distanceUnit; ...
    vInitialPhysical / p.velocityUnit; 1];
targetPlanar = [rFinalPhysical / p.distanceUnit; ...
    vFinalPhysical / p.velocityUnit];

fprintf('Boundary states\n');
fprintf('Initial state [km, km/s]:\n');
fprintf('%+.10f  %+.10f  %+.10f  %+.10f  %+.10f  %+.10f\n', ...
    [rInitialPhysical; vInitialPhysical]);
fprintf('Target state [km, km/s]:\n');
fprintf('%+.10f  %+.10f  %+.10f  %+.10f  %+.10f  %+.10f\n\n', ...
    [rFinalPhysical; vFinalPhysical]);

rhoPhysical = linspace(20, 70, 1200);
qPhysical = dustDensity(rhoPhysical / p.distanceUnit, p) / p.distanceUnit^3;

% Index names by MATLAB figure number so each exported image is meaningful
% when copied outside this script.
figureNames = cell(1, 23);
figureNames(4:15) = {
    'eros_dust_density_profile.png', ...
    'eros_planar_optimal_trajectory.png', ...
    'eros_planar_radius_profile.png', ...
    'eros_planar_thrust_angles.png', ...
    'eros_planar_hamiltonian_history.png', ...
    'eros_inclined_optimal_trajectory.png', ...
    'eros_inclined_radius_profile.png', ...
    'eros_inclined_thrust_angles.png', ...
    'eros_inclined_inclination_history.png', ...
    'eros_inclined_hamiltonian_history.png', ...
    'eros_nbody_trajectory_comparison.png', ...
    'eros_nbody_ntc_error_history.png'};
figureNames(20:23) = {
    'eros_cumulative_dust_exposure.png', ...
    'eros_inclination_trade_space.png', ...
    'eros_perturbation_source_attribution.png', ...
    'eros_osculating_orbit_differences.png'};

figure(4); clf;
set(gcf, 'Position', [100, 100, 538, 420]);
plot(rhoPhysical, qPhysical, 'k-', 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
hold on;
hInitialLine = xline(initialRadius, 'k-.', 'LineWidth', 1.1, ...
    'DisplayName', 'h_i');
hFinalLine = xline(finalRadius, 'k--', 'LineWidth', 1.1, ...
    'DisplayName', 'h_f');
grid on; box on;
xlim([20, 70]);
xlabel('\rho / km');
ylabel('q / 1/km^3');
legend([hInitialLine, hFinalLine], {'h_i', 'h_f'}, 'Location', 'northwest');
set(gca, 'FontSize', 11);

%% Nondimensional mission parameters
fprintf('Nondimensional mission parameters\n');
fprintf('DU = %.10f km, MU = %.10f kg, TU = %.10f s, VU = %.10f km/s\n', ...
    p.distanceUnit, p.massUnit, p.timeUnit, p.velocityUnit);
fprintf('hi = %.10f, hf = %.10f, Ra = %.10f, mu = %.10f\n', ...
    p.hInitial / p.distanceUnit, p.hFinal / p.distanceUnit, ...
    p.asteroidRadius / p.distanceUnit, p.mu);
fprintf('rhoA = %.10f, rhoB = %.10f, m0 = %.10f\n', ...
    p.rhoA, p.rhoB, 1.0);
fprintf('T = %.10f, Isp = %.10f, g0 = %.10f\n\n', ...
    p.thrust, p.isp, p.g0);

%% PMP dynamics and shooting residual
% Library routines in functions/ (canonicalDynamics, canonicalJacobian,
% hamiltonian, shootingResidual) implement the PMP equations and the eight
% terminal conditions [r-rf; v-vf; lambda_m; H] = 0.

%% Planar optimal transfer
numberOfCandidates = 300;
maximumRandomBatches = 5;
% Fixed RNG seed for reproducible costate screening.
randomSeed = 10775298;
[solutionPlanar, timePlanar, canonicalPlanar] = solveRandomTransfer(...
    xInitial, targetPlanar, p, odeOptions, screeningOptions, ...
    numberOfCandidates, maximumRandomBatches, randomSeed);

[positionErrorPlanar, velocityErrorPlanar] = terminalErrors(...
    canonicalPlanar(end, :).', targetPlanar, p);
fprintf('Planar PMP solution\n');
printTransferSolution(solutionPlanar, canonicalPlanar(end, 7), ...
    positionErrorPlanar, velocityErrorPlanar, p);

plotTrajectoryMap(5, timePlanar, canonicalPlanar, p, initialRadius, finalRadius);
plotRadiusProfile(6, timePlanar, canonicalPlanar, p);
plotThrustAngles(7, timePlanar, canonicalPlanar, p);
plotHamiltonian(8, timePlanar, canonicalPlanar, p);

%% Inclination continuation to 3.5 deg
inclinations = linspace(0, 3.5, 5);
solutionInclined = solutionPlanar;
targetInclined = targetPlanar;
inclinationSolutions = cell(size(inclinations));
inclinationTimes = cell(size(inclinations));
inclinationTrajectories = cell(size(inclinations));
inclinationSolutions{1} = solutionPlanar;
inclinationTimes{1} = timePlanar;
inclinationTrajectories{1} = canonicalPlanar;

for index = 2:numel(inclinations)
    angle = deg2rad(inclinations(index));
    targetInclined = [rFinalPhysical / p.distanceUnit; ...
        [0; finalSpeed * cos(angle); finalSpeed * sin(angle)] / p.velocityUnit];
    solutionInclined = continueTransfer(solutionInclined, xInitial, ...
        targetInclined, p, odeOptions);
    [inclinationTimes{index}, inclinationTrajectories{index}] = ...
        propagateCanonical(solutionInclined, xInitial, p, odeOptions, 1800);
    inclinationSolutions{index} = solutionInclined;
    fprintf('Inclination continuation: %.3f deg converged.\n', ...
        inclinations(index));
end
timeInclined = inclinationTimes{end};
canonicalInclined = inclinationTrajectories{end};
[positionErrorInclined, velocityErrorInclined] = terminalErrors(...
    canonicalInclined(end, :).', targetInclined, p);

fprintf('\nInclined PMP solution\n');
printTransferSolution(solutionInclined, canonicalInclined(end, 7), ...
    positionErrorInclined, velocityErrorInclined, p);

plotTrajectoryMap(9, timeInclined, canonicalInclined, p, initialRadius, finalRadius);
plotRadiusProfile(10, timeInclined, canonicalInclined, p);
plotThrustAngles(11, timeInclined, canonicalInclined, p);
plotInclination(12, timeInclined, canonicalInclined, p);
plotHamiltonian(13, timeInclined, canonicalInclined, p);
plotExposureHistory(20, timePlanar, canonicalPlanar, timeInclined, ...
    canonicalInclined, p, initialRadius, finalRadius);
tradeMetrics = plotInclinationTradeSpace(21, inclinations, ...
    inclinationSolutions, inclinationTimes, inclinationTrajectories, p);

% Terminal conditions and Hamiltonian constancy are independent checks.
assert(positionErrorPlanar < 1e-5 && velocityErrorPlanar < 1e-5);
assert(positionErrorInclined < 1e-5 && velocityErrorInclined < 1e-5);

%% SPICE n-body free-flight check
addpath(fullfile(projectDirectory, 'mice', 'src', 'mice'));
addpath(fullfile(projectDirectory, 'mice', 'lib'));
cspice_kclear;
kernelDirectory = fullfile(projectDirectory, 'Kernels');
cspice_furnsh(fullfile(kernelDirectory, 'naif0012.tls'));
cspice_furnsh(fullfile(kernelDirectory, 'de440s.bsp'));
cspice_furnsh(fullfile(kernelDirectory, '2000433.bsp'));
cspice_furnsh(fullfile(kernelDirectory, 'gm_de440.tpc'));
spiceCleanup = onCleanup(@() cspice_kclear);

bodyNames = {'MERCURY BARYCENTER', 'VENUS BARYCENTER', 'EARTH', ...
    'MOON', 'MARS BARYCENTER', 'JUPITER BARYCENTER', ...
    'SATURN BARYCENTER', 'URANUS BARYCENTER', ...
    'NEPTUNE BARYCENTER', 'PLUTO BARYCENTER', 'SUN'};
bodyGm = zeros(size(bodyNames));
for index = 1:numel(bodyNames)
    bodyGm(index) = cspice_bodvrd(bodyNames{index}, 'GM', 1);
end

etInitial = cspice_str2et('2012 JAN 15 00:00:00.000 TDB');
etFinal = cspice_str2et('2012 FEB 15 00:00:00.000 TDB');
propagationTime = linspace(0, etFinal - etInitial, 2500);
finalScaledState = canonicalInclined(end, 1:6).';
finalPhysicalState = [finalScaledState(1:3) * p.distanceUnit; ...
    finalScaledState(4:6) * p.velocityUnit];
nBodyOptions = odeset('RelTol', 1e-12, 'AbsTol', 1e-13);

[timeNBody, stateNBody] = ode113(@(t, x) nBodyDynamics(...
    t, x, etInitial, p.muPhysical, bodyNames, bodyGm), ...
    propagationTime, finalPhysicalState, nBodyOptions);
[~, stateKepler] = ode113(@(t, x) twoBodyDynamics(x, p.muPhysical), ...
    propagationTime, finalPhysicalState, nBodyOptions);

% Isolate the dominant forcing mechanisms. Leave-one-out runs measure each
% body's marginal effect without assuming that nonlinear perturbations add.
sunIndex = find(strcmp(bodyNames, 'SUN'), 1);
stateSun = propagateNBodyHistory(propagationTime, finalPhysicalState, ...
    etInitial, p.muPhysical, bodyNames(sunIndex), bodyGm(sunIndex), nBodyOptions);
nonSunIndices = setdiff(1:numel(bodyNames), sunIndex);
stateNonSun = propagateNBodyHistory(propagationTime, finalPhysicalState, ...
    etInitial, p.muPhysical, bodyNames(nonSunIndices), ...
    bodyGm(nonSunIndices), nBodyOptions);

finalMarginalNtc = zeros(numel(bodyNames), 3);
finalRotation = ntcRotation(stateKepler(end, 1:6).');
for index = 1:numel(bodyNames)
    retained = setdiff(1:numel(bodyNames), index);
    stateWithoutBody = propagateNBodyHistory(propagationTime, ...
        finalPhysicalState, etInitial, p.muPhysical, bodyNames(retained), ...
        bodyGm(retained), nBodyOptions);
    finalMarginalNtc(index, :) = (finalRotation * (...
        stateNBody(end, 1:3) - stateWithoutBody(end, 1:3)).').';
end

plotNBodyTrajectory(14, stateNBody, stateKepler);
plotNtcErrors(15, timeNBody, stateNBody, stateKepler);
plotPerturbationAttribution(22, timeNBody, stateNBody, stateSun, ...
    stateNonSun, stateKepler, bodyNames, finalMarginalNtc);
elementDifferences = plotOsculatingDifferences(23, timeNBody, ...
    stateNBody, stateKepler, p.muPhysical);

fprintf('SPICE n-body free-flight check\n');
fprintf('31-day n-body/Kepler position difference: %.10e km\n', ...
    norm(stateNBody(end, 1:3) - stateKepler(end, 1:3)));
fprintf('31-day n-body/Kepler velocity difference: %.10e km/s\n', ...
    norm(stateNBody(end, 4:6) - stateKepler(end, 4:6)));
fprintf('31-day semimajor-axis difference: %.10e km\n', ...
    elementDifferences(end, 1));
fprintf('31-day phase difference: %.10e deg\n', ...
    elementDifferences(end, 5));
fprintf('Maximum eccentricity difference: %.10e\n', ...
    max(abs(elementDifferences(:, 2))));
fprintf('Planar/inclined nondimensional exposure: %.10e / %.10e\n', ...
    tradeMetrics.exposure(1), tradeMetrics.exposure(end));

assert(all(isfinite(finalMarginalNtc), 'all'));
assert(all(isfinite(elementDifferences), 'all'));
assert(max(abs(elementDifferences(:, 1))) < 1e-2);
assert(max(abs(elementDifferences(:, 2))) < 1e-3);
assert(all(diff(inclinations) > 0));
assert(all(tradeMetrics.finalMass > 0));

% Export the full named figure set produced by this script.
exportNamedFigures(4:15, figureNames, projectDirectory);
exportNamedFigures(20:23, figureNames, projectDirectory);

%% Local functions
function [solution, time, trajectory] = solveRandomTransfer(xInitial, target, ...
        p, odeOptions, screeningOptions, candidateCount, maximumBatches, seed)
% Screen uniform random costates, then solve the best genuine PMP candidates.
rng(seed, 'twister');
initialPeriod = 2 * pi * norm(xInitial(1:3))^(3 / 2);
bestSolution = [];
bestCost = inf;

for batch = 1:maximumBatches
    % Most guesses span the stated interval; a smaller order-one subset
    % resolves the natural Hamiltonian scale. Planar symmetry sets the two
    % out-of-plane costates to zero without prescribing the unknown solution.
    focusedCount = round(candidateCount / 3);
    wideCount = candidateCount - focusedCount;
    costateGuesses = [-100 + 200 * rand(7, wideCount), ...
        -5 + 10 * rand(7, focusedCount)];
    costateGuesses([3, 6], :) = 0;
    guesses = [costateGuesses; ...
        (2 + rand(1, candidateCount)) * initialPeriod];
    scores = inf(1, candidateCount);

    % Cheap integrations rank random guesses before the expensive refinement.
    for index = 1:candidateCount
        try
            residual = shootingResidual(guesses(:, index), xInitial, ...
                target, p, screeningOptions);
            scores(index) = norm(residual);
        catch
            scores(index) = inf;
        end
    end

    [~, ordering] = sort(scores);
    candidatesToSolve = ordering(1:min(12, candidateCount));
    fprintf('Random batch %d: best pre-solve residual %.6e\n', ...
        batch, scores(ordering(1)));

    % Refine only the best candidates from this batch with the full tolerance.
    for candidate = candidatesToSolve
        try
            trial = continueTransfer(guesses(:, candidate), xInitial, ...
                target, p, odeOptions);
            [trialTime, trialTrajectory] = propagateCanonical(...
                trial, xInitial, p, odeOptions, 1200);
            trialResidual = shootingResidual(trial, xInitial, target, ...
                p, odeOptions);
            trialCost = trapz(trialTime, dustDensity(...
                vecnorm(trialTrajectory(:, 1:3), 2, 2), p));

            if norm(trialResidual) < 1e-8 && trialCost < bestCost
                bestCost = trialCost;
                bestSolution = trial;
            end
        catch
            % A failed random candidate is discarded; the next candidate is independent.
        end
    end

    if ~isempty(bestSolution)
        solution = bestSolution;
        [time, trajectory] = propagateCanonical(...
            solution, xInitial, p, odeOptions, 1800);
        return;
    end
end

error('No converged PMP extremal was found after %d random batches.', ...
    maximumBatches);
end

function solution = continueTransfer(initialGuess, xInitial, target, p, odeOptions)
% Solve the square eight-condition shooting system with analytical derivatives.
solverOptions = optimoptions('fsolve', ...
    'Algorithm', 'trust-region-dogleg', ...
    'SpecifyObjectiveGradient', true, ...
    'Display', 'off', ...
    'FunctionTolerance', 1e-11, ...
    'StepTolerance', 1e-12, ...
    'OptimalityTolerance', 1e-11, ...
    'MaxIterations', 250, ...
    'MaxFunctionEvaluations', 3000);

[solution, residual, exitFlag] = fsolve(@(z) shootingResidual(...
    z, xInitial, target, p, odeOptions), initialGuess, solverOptions);

if exitFlag <= 0 || norm(residual) > 1e-8 || solution(8) <= 0
    error('Single shooting failed: exit flag %d, residual %.3e.', ...
        exitFlag, norm(residual));
end
end

function [residual, jacobian] = shootingResidual(decision, xInitial, ...
        target, p, odeOptions)
% Propagate the canonical system and evaluate terminal constraints and Jacobian.
canonicalInitial = [xInitial; decision(1:7)];

if nargout < 2
    [~, trajectory] = ode113(@(t, y) canonicalDynamics(y, p), ...
        [0, decision(8)], canonicalInitial, odeOptions);
    residual = boundaryResidual(trajectory(end, :).', target, p);
    return;
end

sensitivityInitial = [zeros(7, 7); eye(7)];
augmentedInitial = [canonicalInitial; sensitivityInitial(:)];
[~, augmented] = ode113(@(t, y) canonicalVariational(y, p), ...
    [0, decision(8)], augmentedInitial, odeOptions);
finalCanonical = augmented(end, 1:14).';
finalSensitivity = reshape(augmented(end, 15:end), 14, 7);
[residual, residualStateJacobian] = boundaryResidual(...
    finalCanonical, target, p);
finalDerivative = canonicalDynamics(finalCanonical, p);
jacobian = [residualStateJacobian * finalSensitivity, ...
    residualStateJacobian * finalDerivative];
end

function [residual, stateJacobian] = boundaryResidual(canonical, target, p)
% Enforce final position, velocity, free mass, and free-time transversality.
residual = [canonical(1:6) - target; canonical(14); ...
    hamiltonian(canonical, p)];

if nargout > 1
    stateJacobian = zeros(8, 14);
    stateJacobian(1:6, 1:6) = eye(6);
    stateJacobian(7, 14) = 1;
    stateJacobian(8, :) = hamiltonianGradient(canonical, p).';
end
end

function [time, trajectory] = propagateCanonical(decision, xInitial, p, ...
        odeOptions, sampleCount)
% Reintegrate a converged extremal on a uniform grid for diagnostics and plots.
timeGrid = linspace(0, decision(8), sampleCount);
[time, trajectory] = ode113(@(t, y) canonicalDynamics(y, p), ...
    timeGrid, [xInitial; decision(1:7)], odeOptions);
end

function derivative = canonicalVariational(augmented, p)
% Propagate the canonical state and its sensitivity to the initial costates.
canonical = augmented(1:14);
sensitivity = reshape(augmented(15:end), 14, 7);
derivative = [canonicalDynamics(canonical, p); ...
    reshape(canonicalJacobian(canonical, p) * sensitivity, 98, 1)];
end

function derivative = canonicalDynamics(canonical, p)
% Evaluate the full-thrust PMP state and costate equations.
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

function jacobian = canonicalJacobian(canonical, p)
% Evaluate the analytical Jacobian used by the shooting sensitivities.
r = canonical(1:3);
m = canonical(7);
lambdaV = canonical(11:13);
radius = norm(r);
primerNorm = norm(lambdaV);
alphaDirection = lambdaV / primerNorm;
gravityMatrix = eye(3) / radius^3 - 3 * (r * r.') / radius^5;
primerDerivative = eye(3) / primerNorm - ...
    (lambdaV * lambdaV.') / primerNorm^3;

scalarProduct = dot(r, lambdaV);
gravityCostateJacobian = -3 * (...
    lambdaV * r.' + r * lambdaV.' + scalarProduct * eye(3)) / radius^5 ...
    + 15 * scalarProduct * (r * r.') / radius^7;
[~, ~, dustHessian] = dustDerivatives(radius, r, p);

jacobian = zeros(14);
jacobian(1:3, 4:6) = eye(3);
jacobian(4:6, 1:3) = -gravityMatrix;
jacobian(4:6, 7) = p.thrust * alphaDirection / m^2;
jacobian(4:6, 11:13) = -p.thrust * primerDerivative / m;
jacobian(8:10, 1:3) = gravityCostateJacobian - dustHessian;
jacobian(8:10, 11:13) = gravityMatrix;
jacobian(11:13, 8:10) = -eye(3);
jacobian(14, 7) = 2 * p.thrust * primerNorm / m^3;
jacobian(14, 11:13) = -p.thrust * alphaDirection.' / m^2;
end

function value = hamiltonian(canonical, p)
% Evaluate the autonomous minimum-principle Hamiltonian.
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

function gradient = hamiltonianGradient(canonical, p)
% Evaluate dH/d[r,v,m,lambda_r,lambda_v,lambda_m] analytically.
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

function q = dustDensity(radius, p)
% Evaluate nondimensional radial dust density for scalar or vector radii.
q = p.k1 ./ (p.k2 + (radius - p.rhoA).^2) ...
    + p.k3 ./ (p.k4 + (radius - p.rhoB).^2);
end

function [q, gradient, hessian] = dustDerivatives(radius, r, p)
% Evaluate a radial density together with its Cartesian gradient and Hessian.
dA = radius - p.rhoA;
dB = radius - p.rhoB;
denominatorA = p.k2 + dA^2;
denominatorB = p.k4 + dB^2;
q = p.k1 / denominatorA + p.k3 / denominatorB;
qPrime = -2 * p.k1 * dA / denominatorA^2 ...
    - 2 * p.k3 * dB / denominatorB^2;
radialDirection = r / radius;
gradient = qPrime * radialDirection;

if nargout > 2
    qSecond = -2 * p.k1 / denominatorA^2 ...
        + 8 * p.k1 * dA^2 / denominatorA^3 ...
        - 2 * p.k3 / denominatorB^2 ...
        + 8 * p.k3 * dB^2 / denominatorB^3;
    hessian = qSecond * (radialDirection * radialDirection.') ...
        + (qPrime / radius) * (eye(3) ...
        - radialDirection * radialDirection.');
end
end

function [positionError, velocityError] = terminalErrors(finalCanonical, target, p)
% Convert terminal residuals to km and m/s.
positionError = norm(finalCanonical(1:3) - target(1:3)) * p.distanceUnit;
velocityError = norm(finalCanonical(4:6) - target(4:6)) ...
    * p.velocityUnit * 1000;
end

function printTransferSolution(solution, finalMass, positionError, velocityError, p)
% Print the tabulated transfer quantities at high precision.
fprintf('tf = %.10f min, mf = %.10f kg\n', ...
    solution(8) * p.timeUnit / 60, finalMass * p.massUnit);
fprintf('lambda0 =');
fprintf(' %+.10f', solution(1:7));
fprintf('\nPosition error = %.10e km\n', positionError);
fprintf('Velocity error = %.10e m/s\n\n', velocityError);
end

function plotTrajectoryMap(number, ~, canonical, p, initialRadius, finalRadius)
% Overlay the optimal trajectory on the physical dust-density field.
limits = linspace(-80, 80, 400);
[xGrid, yGrid] = meshgrid(limits, limits);
radiusGrid = hypot(xGrid, yGrid);
densityGrid = dustDensity(radiusGrid / p.distanceUnit, p) / p.distanceUnit^3;
position = canonical(:, 1:3) * p.distanceUnit;
angle = linspace(0, 2 * pi, 600);

figure(number); clf;
if number == 5
    set(gcf, 'Position', [100, 100, 686, 564]);
else
    set(gcf, 'Position', [100, 100, 654, 564]);
end
surface(xGrid, yGrid, zeros(size(xGrid)), densityGrid, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
view(2); hold on;
plot(position(:, 1), position(:, 2), 'k-', 'LineWidth', 1.5);
plot(initialRadius * cos(angle), initialRadius * sin(angle), ...
    'k:', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(finalRadius * cos(angle), finalRadius * sin(angle), ...
    'k--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
startHandle = plot(position(1, 1), position(1, 2), 'ko', 'MarkerSize', 9, ...
    'LineWidth', 1.2, 'DisplayName', 'Start');
endHandle = plot(position(end, 1), position(end, 2), 'kx', 'MarkerSize', 9, ...
    'LineWidth', 1.2, 'DisplayName', 'End');
axis equal; xlim([-80, 80]); ylim([-80, 80]);
grid on; box on; colormap(parula);
colorbarHandle = colorbar;
colorbarHandle.Label.String = 'q / 1/km^3';
xlabel('x / km'); ylabel('y / km');
legend([startHandle, endHandle], {'Start', 'End'}, 'Location', 'northwest');
set(gca, 'FontSize', 11);
end

function plotRadiusProfile(number, time, canonical, p)
% Compare the radial history with the outer density peak.
radiusPhysical = vecnorm(canonical(:, 1:3), 2, 2) * p.distanceUnit;
rho = linspace(52, 70, 800);
densityMicro = dustDensity(rho / p.distanceUnit, p) ...
    / p.distanceUnit^3 * 1e6;

figure(number); clf;
if number == 6
    set(gcf, 'Position', [100, 100, 956, 428]);
else
    set(gcf, 'Position', [100, 100, 944, 432]);
end
layout = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout, 1);
plot(densityMicro, rho, 'k-', 'LineWidth', 1.4); hold on;
yline(p.rhoBPhysical, 'k--', 'LineWidth', 1.0);
xlim([0, 1.05 * max(densityMicro)]); ylim([52, 70]);
grid on; box on;
xlabel('q / (10^{-6} km^{-3})'); ylabel('\rho / km');

nexttile(layout, 2, [1, 3]);
plot(time * p.timeUnit, radiusPhysical, 'k-', 'LineWidth', 1.4); hold on;
yline(p.rhoBPhysical, 'k--', 'LineWidth', 1.0);
xlim([0, time(end) * p.timeUnit]); ylim([52, 70]);
grid on; box on;
xlabel('t / s'); ylabel('\rho / km');
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
end

function plotThrustAngles(number, time, canonical, p)
% Plot thrust direction in the instantaneous radial-transverse-cross frame.
[inPlane, outOfPlane] = thrustAngles(canonical);
figure(number); clf;
if number == 7
    set(gcf, 'Position', [100, 100, 956, 570]);
else
    set(gcf, 'Position', [100, 100, 952, 566]);
end
layout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout);
plot(time * p.timeUnit, inPlane, 'k-', 'LineWidth', 1.4);
grid on; box on;
ylabel('\alpha_{T,plane} / deg');
xlim([0, time(end) * p.timeUnit]);
nexttile(layout);
plot(time * p.timeUnit, outOfPlane, 'k-', 'LineWidth', 1.4);
grid on; box on;
xlabel('t / s'); ylabel('\alpha_{T,cross} / deg');
xlim([0, time(end) * p.timeUnit]);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
end

function [inPlane, outOfPlane] = thrustAngles(canonical)
% Resolve the primer control relative to retrograde, radial, and orbit-normal axes.
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

function plotInclination(number, time, canonical, p)
% Compute osculating inclination from the angular-momentum direction.
angularMomentum = cross(canonical(:, 1:3), canonical(:, 4:6), 2);
inclination = acosd(angularMomentum(:, 3) ...
    ./ vecnorm(angularMomentum, 2, 2));
figure(number); clf;
set(gcf, 'Position', [100, 100, 594, 390]);
plot(time * p.timeUnit, inclination, 'k-', 'LineWidth', 1.4);
grid on; box on;
xlabel('t / s'); ylabel('inc / deg');
xlim([0, time(end) * p.timeUnit]); ylim([0, 3.5]);
set(gca, 'FontSize', 11);
end

function plotHamiltonian(number, time, canonical, p)
% Confirm the autonomous free-time Hamiltonian remains constant and near zero.
values = zeros(size(time));
for index = 1:numel(time)
    values(index) = hamiltonian(canonical(index, :).', p);
end
figure(number); clf;
if number == 8
    set(gcf, 'Position', [100, 100, 600, 436]);
else
    set(gcf, 'Position', [100, 100, 552, 446]);
end
plot(time * p.timeUnit, values, 'k-', 'LineWidth', 1.4);
grid on; box on;
xlabel('t / s'); ylabel('Hamiltonian');
xlim([0, time(end) * p.timeUnit]);
set(gca, 'FontSize', 11);
end

function plotExposureHistory(number, timePlanar, planar, timeInclined, ...
        inclined, p, initialRadius, finalRadius)
% Show instantaneous and accumulated dust exposure along both extremals.
radiusPlanar = vecnorm(planar(:, 1:3), 2, 2);
radiusInclined = vecnorm(inclined(:, 1:3), 2, 2);
densityPlanar = dustDensity(radiusPlanar, p);
densityInclined = dustDensity(radiusInclined, p);
cumulativePlanar = cumtrapz(timePlanar, densityPlanar);
cumulativeInclined = cumtrapz(timeInclined, densityInclined);

% This same-duration monotonic-radius path is explicitly kinematic: it is an
% exposure yardstick, not a dynamically feasible replacement trajectory.
referenceRadius = linspace(initialRadius / p.distanceUnit, ...
    finalRadius / p.distanceUnit, numel(timePlanar)).';
referenceDensity = dustDensity(referenceRadius, p);
referenceCumulative = cumtrapz(timePlanar, referenceDensity);

figure(number); clf;
set(gcf, 'Position', [100, 100, 850, 620]);
layout = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout);
plot(timePlanar / timePlanar(end), densityPlanar, 'k-', 'LineWidth', 1.3);
hold on;
plot(timeInclined / timeInclined(end), densityInclined, 'Color', ...
    [0.85, 0.325, 0.098], 'LineWidth', 1.2);
plot(timePlanar / timePlanar(end), referenceDensity, 'k--', 'LineWidth', 1.0);
grid on; box on; ylabel('q(\rho) / ND');
legend({'Planar optimum', 'Inclined optimum', ...
    'Monotonic-radius yardstick'}, 'Location', 'northwest');
nexttile(layout);
plot(timePlanar / timePlanar(end), cumulativePlanar, 'k-', 'LineWidth', 1.3);
hold on;
plot(timeInclined / timeInclined(end), cumulativeInclined, 'Color', ...
    [0.85, 0.325, 0.098], 'LineWidth', 1.2);
plot(timePlanar / timePlanar(end), referenceCumulative, 'k--', 'LineWidth', 1.0);
grid on; box on; xlabel('t / t_f'); ylabel('\int_0^t q(\rho) d\tau / ND');
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
end

function metrics = plotInclinationTradeSpace(number, inclinations, solutions, ...
        times, trajectories, p)
% Quantify how exposure, duration, propellant, and steering vary with plane change.
sampleCount = numel(inclinations);
metrics = struct('exposure', zeros(1, sampleCount), ...
    'tofMinutes', zeros(1, sampleCount), ...
    'finalMass', zeros(1, sampleCount), ...
    'peakDensity', zeros(1, sampleCount), ...
    'meanCrossThrust', zeros(1, sampleCount));

for index = 1:sampleCount
    time = times{index};
    canonical = trajectories{index};
    density = dustDensity(vecnorm(canonical(:, 1:3), 2, 2), p);
    primer = canonical(:, 11:13);
    thrustDirection = -primer ./ vecnorm(primer, 2, 2);
    metrics.exposure(index) = trapz(time, density);
    metrics.tofMinutes(index) = solutions{index}(8) * p.timeUnit / 60;
    metrics.finalMass(index) = canonical(end, 7) * p.massUnit;
    metrics.peakDensity(index) = max(density);
    metrics.meanCrossThrust(index) = trapz(time, ...
        abs(thrustDirection(:, 3))) / time(end);
end

figure(number); clf;
set(gcf, 'Position', [100, 100, 900, 700]);
layout = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
plotTradeTile(layout, inclinations, metrics.exposure, ...
    'J = \int q d\tau / ND');
plotTradeTile(layout, inclinations, metrics.tofMinutes, 't_f / min');
plotTradeTile(layout, inclinations, p.initialMassPhysical - metrics.finalMass, ...
    'Propellant / kg');
plotTradeTile(layout, inclinations, metrics.peakDensity, 'max q / ND');
plotTradeTile(layout, inclinations, metrics.meanCrossThrust, ...
    'mean |\alpha_C|');
nexttile(layout); axis off;
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
end

function plotTradeTile(layout, inclination, value, label)
% Draw one consistently formatted inclination trade-space quantity.
nexttile(layout);
plot(inclination, value, 'ko-', 'MarkerFaceColor', 'w', 'LineWidth', 1.2);
grid on; box on; xlabel('\Delta i / deg'); ylabel(label);
end

function state = propagateNBodyHistory(time, initialState, etInitial, ...
        muEros, names, gm, odeOptions)
% Propagate a selected subset of third bodies on the common output grid.
[~, state] = ode113(@(t, x) nBodyDynamics(t, x, etInitial, ...
    muEros, names, gm), time, initialState, odeOptions);
end

function derivative = nBodyDynamics(time, state, etInitial, muEros, names, gm)
% Evaluate Eros gravity plus differential third-body accelerations from SPICE.
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

function derivative = twoBodyDynamics(state, mu)
% Evaluate the unperturbed Eros-centered Kepler reference dynamics.
derivative = [state(4:6); -mu * state(1:3) / norm(state(1:3))^3];
end

function plotNBodyTrajectory(number, nBody, kepler)
% Show the ecliptic-plane and side views of the perturbed inclined orbit.
figure(number); clf;
set(gcf, 'Position', [100, 100, 592, 480]);
layout = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout, 1, [3, 1]);
plot(kepler(:, 1), kepler(:, 2), 'Color', [0.65, 0.65, 0.65], ...
    'LineStyle', '--', 'LineWidth', 1.0); hold on;
plot(nBody(:, 1), nBody(:, 2), 'k-', 'LineWidth', 1.2);
axis equal; grid on; box on;
xlabel('x / km'); ylabel('y / km');
xlim([-60, 60]); ylim([-60, 60]);
nexttile(layout, 4);
plot(nBody(:, 2), nBody(:, 3), 'k-', 'LineWidth', 1.2);
grid on; box on;
xlabel('y / km'); ylabel('z / km');
xlim([-60, 60]);
ylim([-5, 5]);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 10);
end

function plotNtcErrors(number, time, nBody, kepler)
% Resolve n-body minus Kepler errors in radial (N), tangential (T), and cross (C).
[positionNtc, velocityNtc, positionDifference, velocityDifference] = ...
    ntcDifferences(nBody, kepler);

days = time / 86400;
figure(number); clf;
set(gcf, 'Position', [100, 100, 966, 698]);
layout = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
positionLabels = {'|\Delta r_N| / km', '|\Delta r_T| / km', ...
    '|\Delta r_C| / km'};
velocityLabels = {'|\Delta v_N| / km/s', '|\Delta v_T| / km/s', ...
    '|\Delta v_C| / km/s'};

for component = 1:3
    nexttile(layout, 2 * component - 1);
    semilogy(days, max(abs(positionNtc(:, component)), 1e-12), ...
        'k-', 'LineWidth', 1.0);
    grid on; box on; ylabel(positionLabels{component}); xlim([0, days(end)]);
    nexttile(layout, 2 * component);
    semilogy(days, max(abs(velocityNtc(:, component)), 1e-12), ...
        'k-', 'LineWidth', 1.0);
    grid on; box on; ylabel(velocityLabels{component}); xlim([0, days(end)]);
end

nexttile(layout, 7);
semilogy(days, max(vecnorm(positionDifference, 2, 2), 1e-12), ...
    'k-', 'LineWidth', 1.0);
grid on; box on; xlabel('t / day'); ylabel('||\Delta r_{in}|| / km');
xlim([0, days(end)]);
nexttile(layout, 8);
semilogy(days, max(vecnorm(velocityDifference, 2, 2), 1e-12), ...
    'k-', 'LineWidth', 1.0);
grid on; box on; xlabel('t / day'); ylabel('||\Delta v_{in}|| / km/s');
xlim([0, days(end)]);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 8);
end

function plotPerturbationAttribution(number, time, fullState, sunState, ...
        nonSunState, kepler, bodyNames, finalMarginalNtc)
% Compare signed errors and each body's leave-one-out marginal contribution.
[positionFull, velocityFull] = ntcDifferences(fullState, kepler);
[positionSun, ~] = ntcDifferences(sunState, kepler);
[positionNonSun, ~] = ntcDifferences(nonSunState, kepler);
days = time / 86400;

figure(number); clf;
set(gcf, 'Position', [100, 100, 1000, 720]);
layout = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(layout);
plot(days, positionFull, 'LineWidth', 1.0);
grid on; box on; xlabel('t / day'); ylabel('\Delta r_{NTC} / km');
legend({'N', 'T', 'C'}, 'Location', 'northwest');
nexttile(layout);
plot(days, velocityFull * 1000, 'LineWidth', 1.0);
grid on; box on; xlabel('t / day'); ylabel('\Delta v_{NTC} / m/s');
legend({'N', 'T', 'C'}, 'Location', 'northwest');
nexttile(layout);
plot(days, positionFull(:, 2), 'k-', 'LineWidth', 1.3); hold on;
plot(days, positionSun(:, 2), '--', 'Color', [0.85, 0.325, 0.098], ...
    'LineWidth', 1.2);
plot(days, positionNonSun(:, 2), 'k:', 'LineWidth', 1.2);
grid on; box on; xlabel('t / day'); ylabel('\Delta r_T / km');
legend({'All bodies', 'Sun only', 'All except Sun'}, 'Location', 'northwest');
nexttile(layout);
bar(finalMarginalNtc(:, 2), 'FaceColor', [0.35, 0.35, 0.35]);
grid on; box on; ylabel('Final marginal \Delta r_T / km');
labels = strrep(bodyNames, ' BARYCENTER', '');
set(gca, 'XTick', 1:numel(labels), 'XTickLabel', labels);
xtickangle(45);
set(findall(gcf, 'Type', 'axes'), 'FontSize', 9);
end

function differences = plotOsculatingDifferences(number, time, nBody, ...
        kepler, mu)
% Use nonsingular phase diagnostics for the nearly circular inclined orbit.
nBodyElements = orbitalDiagnostics(nBody, mu);
keplerElements = orbitalDiagnostics(kepler, mu);
differences = nBodyElements - keplerElements;
differences(:, 3:5) = rad2deg(differences(:, 3:5));
days = time / 86400;

% The two-body reference must preserve its energy-derived semimajor axis.
assert(max(abs(keplerElements(:, 1) - keplerElements(1, 1))) < 1e-6);

figure(number); clf;
set(gcf, 'Position', [100, 100, 900, 700]);
layout = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
labels = {'\Delta a / km', '\Delta e', '\Delta i / deg', ...
    '\Delta \Omega / deg', '\Delta u / deg'};
for index = 1:5
    nexttile(layout);
    plot(days, differences(:, index), 'k-', 'LineWidth', 1.1);
    grid on; box on; xlabel('t / day'); ylabel(labels{index});
end
nexttile(layout); axis off;
set(findall(gcf, 'Type', 'axes'), 'FontSize', 9);
end

function elements = orbitalDiagnostics(state, mu)
% Return [a,e,i,Omega,u]; argument of latitude u avoids singular periapsis angles.
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

function [positionNtc, velocityNtc, positionDifference, velocityDifference] = ...
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

function rotation = ntcRotation(state)
% Construct the radial-tangential-cross rotation from an inertial state.
normal = state(1:3) / norm(state(1:3));
crossTrack = cross(state(1:3), state(4:6));
crossTrack = crossTrack / norm(crossTrack);
tangential = cross(crossTrack, normal);
rotation = [normal, tangential, crossTrack].';
end

function exportNamedFigures(numbers, figureNames, projectDirectory)
% Export named figures to plots/, or to SGN_OUTPUT_DIR when that env var is set.
outputDirectory = getenv('SGN_OUTPUT_DIR');
if isempty(outputDirectory)
    outputDirectory = fullfile(projectDirectory, 'plots');
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
for number = numbers
    % Fail early if a new figure has no meaningful external file name.
    assert(number <= numel(figureNames) && ~isempty(figureNames{number}), ...
        'Missing export file name for figure %d.', number);
    drawnow;
    exportgraphics(figure(number), fullfile(outputDirectory, ...
        figureNames{number}), 'Resolution', 160);
end
end
