% eros_continuous_guidance.m
% Eros fixed-throttle low-thrust transfer: PMP single shooting,
% inclination continuation, and SPICE third-body sensitivity analysis.
% Author: Pasquale Marzaioli

clear; close all; clc;
projectDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDirectory, 'functions'));

%% Shared physical data and nondimensionalization
p.hInitial = 52.25;          % km
p.hFinal = 36.15;            % km
p.referenceRadiusPhysical = 17.00; % km, chosen spherical radius convention
% Match the Eros GM distributed in the bundled DE440 constants kernel.
p.muPhysical = 4.463e-4;     % km^3/s^2
initializationMuPhysical = 3.5e-4;
initializationThrustPhysical = 21.59e-9; % kg km/s^2
p.rhoAPhysical = 23.314 + p.referenceRadiusPhysical;
p.rhoBPhysical = 42.170 + p.referenceRadiusPhysical;
p.k1 = 7.393750e-3;
p.k2 = 7.500000e-3;
p.k3 = 3.696875e-4;
p.k4 = 6.250000e-4;
p.initialMassPhysical = 25.0; % kg
% Preserve the original nondimensional thrust-to-gravity ratio after the GM
% correction; the resulting 27.53 micro-newton thruster is a study assumption.
p.thrustPhysical = initializationThrustPhysical ...
    * p.muPhysical / initializationMuPhysical;
p.ispPhysical = 382.82;       % s
p.g0Physical = 9.80665e-3;    % km/s^2
p.distanceUnit = p.hInitial + p.referenceRadiusPhysical;
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

%% Synthetic radial cost and circular boundary states
initialRadius = p.hInitial + p.referenceRadiusPhysical;
finalRadius = p.hFinal + p.referenceRadiusPhysical;
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
qProxy = dustDensity(rhoPhysical / p.distanceUnit, p);

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
plot(rhoPhysical, qProxy, 'k-', 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');
hold on;
hInitialLine = xline(initialRadius, 'k-.', 'LineWidth', 1.1, ...
    'DisplayName', 'r_i');
hFinalLine = xline(finalRadius, 'k--', 'LineWidth', 1.1, ...
    'DisplayName', 'r_f');
grid on; box on;
xlim([20, 70]);
xlabel('\rho / km');
ylabel('q(\rho)');
legend([hInitialLine, hFinalLine], {'r_i', 'r_f'}, 'Location', 'northwest');
set(gca, 'FontSize', 11);

%% Nondimensional mission parameters
fprintf('Physical mission parameters\n');
fprintf(['Eros GM = %.10e km^3/s^2, reference radius = %.4f km, ' ...
    'thrust = %.6f micro-newton\n'], p.muPhysical, ...
    p.referenceRadiusPhysical, p.thrustPhysical * 1e9);
fprintf('Initial mass = %.6f kg, Isp = %.6f s\n\n', ...
    p.initialMassPhysical, p.ispPhysical);

fprintf('Nondimensional mission parameters\n');
fprintf('DU = %.10f km, MU = %.10f kg, TU = %.10f s, VU = %.10f km/s\n', ...
    p.distanceUnit, p.massUnit, p.timeUnit, p.velocityUnit);
fprintf('hi = %.10f, hf = %.10f, Ra = %.10f, mu = %.10f\n', ...
    p.hInitial / p.distanceUnit, p.hFinal / p.distanceUnit, ...
    p.referenceRadiusPhysical / p.distanceUnit, p.mu);
fprintf('rhoA = %.10f, rhoB = %.10f, m0 = %.10f\n', ...
    p.rhoA, p.rhoB, 1.0);
fprintf('T = %.10f, Isp = %.10f, g0 = %.10f\n\n', ...
    p.thrust, p.isp, p.g0);

%% PMP dynamics and shooting residual
% Library routines in functions/ (canonicalDynamics, canonicalJacobian,
% hamiltonian, shootingResidual) implement the PMP equations and the eight
% terminal conditions [r-rf; v-vf; lambda_m; H] = 0.

%% Planar fixed-throttle extremal
numberOfCandidates = 300;
maximumRandomBatches = 2;
% Fixed RNG seed for reproducible costate screening.
randomSeed = 10775298;

% The measured-GM problem has a narrow shooting basin. Find a reproducible
% anchor with the documented legacy scaling, then continue the coupled physical
% GM/thrust scaling to the corrected model. The legacy model is not a result.
initializationModel = p;
initializationModel.muPhysical = initializationMuPhysical;
initializationModel.thrustPhysical = initializationThrustPhysical;
initializationModel.timeUnit = sqrt(p.distanceUnit^3 ...
    / initializationModel.muPhysical);
initializationModel.velocityUnit = p.distanceUnit ...
    / initializationModel.timeUnit;
initializationModel.thrust = initializationModel.thrustPhysical / (p.massUnit ...
    * p.distanceUnit / initializationModel.timeUnit^2);
initializationModel.isp = p.ispPhysical / initializationModel.timeUnit;
initializationModel.g0 = p.g0Physical * initializationModel.timeUnit^2 ...
    / p.distanceUnit;
initializationModel.massRate = -initializationModel.thrust ...
    / (initializationModel.isp * initializationModel.g0);

[solutionPlanar, ~, ~] = solveRandomTransfer(...
    xInitial, targetPlanar, initializationModel, odeOptions, screeningOptions, ...
    numberOfCandidates, maximumRandomBatches, randomSeed);

muContinuation = linspace(initializationModel.muPhysical, p.muPhysical, 33);
for index = 2:numel(muContinuation)
    continuationModel = p;
    continuationModel.muPhysical = muContinuation(index);
    continuationModel.thrustPhysical = initializationThrustPhysical ...
        * continuationModel.muPhysical / initializationMuPhysical;
    continuationModel.timeUnit = sqrt(p.distanceUnit^3 ...
        / continuationModel.muPhysical);
    continuationModel.velocityUnit = p.distanceUnit ...
        / continuationModel.timeUnit;
    continuationModel.thrust = continuationModel.thrustPhysical / (p.massUnit ...
        * p.distanceUnit / continuationModel.timeUnit^2);
    continuationModel.isp = p.ispPhysical / continuationModel.timeUnit;
    continuationModel.g0 = p.g0Physical * continuationModel.timeUnit^2 ...
        / p.distanceUnit;
    continuationModel.massRate = -continuationModel.thrust ...
        / (continuationModel.isp * continuationModel.g0);
    solutionPlanar = continueTransfer(solutionPlanar, xInitial, targetPlanar, ...
        continuationModel, odeOptions);
    fprintf('GM continuation: %.6e km^3/s^2 converged.\n', ...
        continuationModel.muPhysical);
end
[timePlanar, canonicalPlanar] = propagateCanonical(...
    solutionPlanar, xInitial, p, odeOptions, 1800);

[positionErrorPlanar, velocityErrorPlanar] = terminalErrors(...
    canonicalPlanar(end, :).', targetPlanar, p);
fprintf('Planar PMP solution\n');
printTransferSolution(solutionPlanar, canonicalPlanar(end, 7), ...
    positionErrorPlanar, velocityErrorPlanar, p);

plotTrajectoryMap(5, timePlanar, canonicalPlanar, p, initialRadius, finalRadius);
plotRadiusProfile(6, timePlanar, canonicalPlanar, p);
plotThrustAngles(7, timePlanar, canonicalPlanar, p);
hamiltonianPlanar = plotHamiltonian(8, timePlanar, canonicalPlanar, p);

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
hamiltonianInclined = plotHamiltonian(13, timeInclined, canonicalInclined, p);
plotExposureHistory(20, timePlanar, canonicalPlanar, timeInclined, ...
    canonicalInclined, p, initialRadius, finalRadius);
tradeMetrics = plotInclinationTradeSpace(21, inclinations, ...
    inclinationSolutions, inclinationTimes, inclinationTrajectories, p);

% Terminal conditions and Hamiltonian constancy are independent checks.
assert(positionErrorPlanar < 1e-5 && velocityErrorPlanar < 1e-5);
assert(positionErrorInclined < 1e-5 && velocityErrorInclined < 1e-5);
assert(max(abs(hamiltonianPlanar)) < 1e-8);
assert(max(abs(hamiltonianInclined)) < 1e-8);
fprintf('Maximum |H|, planar/inclined: %.10e / %.10e\n\n', ...
    max(abs(hamiltonianPlanar)), max(abs(hamiltonianInclined)));

%% SPICE third-body free-flight sensitivity check
addpath(fullfile(projectDirectory, 'mice', 'src', 'mice'));
addpath(fullfile(projectDirectory, 'mice', 'lib'));
cspice_kclear;
kernelDirectory = fullfile(projectDirectory, 'Kernels');
cspice_furnsh(fullfile(kernelDirectory, 'naif0012.tls'));
cspice_furnsh(fullfile(kernelDirectory, 'de440s.bsp'));
cspice_furnsh(fullfile(kernelDirectory, '2000433.bsp'));
cspice_furnsh(fullfile(kernelDirectory, 'gm_de440.tpc'));
spiceCleanup = onCleanup(@() cspice_kclear);

% Prevent the guidance and ephemeris checks from silently using different GM values.
muErosKernel = cspice_bodvrd('EROS', 'GM', 1);
assert(abs(muErosKernel - p.muPhysical) / muErosKernel < 1e-12);

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

% Compare the dominant Sun term with the aggregate non-Sun term without eleven
% expensive leave-one-out propagations.
sunIndex = find(strcmp(bodyNames, 'SUN'), 1);
stateSun = propagateNBodyHistory(propagationTime, finalPhysicalState, ...
    etInitial, p.muPhysical, bodyNames(sunIndex), bodyGm(sunIndex), nBodyOptions);
nonSunIndices = setdiff(1:numel(bodyNames), sunIndex);
stateNonSun = propagateNBodyHistory(propagationTime, finalPhysicalState, ...
    etInitial, p.muPhysical, bodyNames(nonSunIndices), ...
    bodyGm(nonSunIndices), nBodyOptions);

plotNBodyTrajectory(14, stateNBody, stateKepler);
plotNtcErrors(15, timeNBody, stateNBody, stateKepler);
plotPerturbationAttribution(22, timeNBody, stateNBody, stateSun, ...
    stateNonSun, stateKepler);
elementDifferences = plotOsculatingDifferences(23, timeNBody, ...
    stateNBody, stateKepler, p.muPhysical);

fprintf('SPICE third-body sensitivity check\n');
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

assert(all(isfinite(elementDifferences), 'all'));
assert(max(abs(elementDifferences(:, 1))) < 1e-2);
assert(max(abs(elementDifferences(:, 2))) < 1e-3);
assert(all(diff(inclinations) > 0));
assert(all(tradeMetrics.finalMass > 0));

% Export the full named figure set produced by this script.
exportNamedFigures(4:15, figureNames, projectDirectory);
exportNamedFigures(20:23, figureNames, projectDirectory);
