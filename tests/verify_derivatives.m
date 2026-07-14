% verify_derivatives.m — Regression checks for analytic derivatives and NTC transforms.
% Uses central finite differences and identity comparisons without extra test frameworks.
% Author: Pasquale Marzaioli

clear;
projectDirectory = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectDirectory, 'functions'));

% Use a nonsingular canonical state representative of the solved trajectories.
p.rhoA = 40.314 / 69.25;
p.rhoB = 59.170 / 69.25;
p.k1 = 7.393750e-3;
p.k2 = 7.500000e-3;
p.k3 = 3.696875e-4;
p.k4 = 6.250000e-4;
p.thrust = 0.0118327079;
p.massRate = -p.thrust / (0.0140339006 * 105373.9704024759);
canonical = [0.74; -0.19; 0.11; 0.03; 0.92; -0.22; 0.91; ...
    0.41; -0.33; 0.27; -0.52; 0.71; 0.38; 0.29];

% Verify the canonical Jacobian column by column.
analyticJacobian = canonicalJacobian(canonical, p);
finiteJacobian = zeros(14);
for column = 1:14
    step = 1e-6 * max(1, abs(canonical(column)));
    plus = canonical;
    minus = canonical;
    plus(column) = plus(column) + step;
    minus(column) = minus(column) - step;
    finiteJacobian(:, column) = (canonicalDynamics(plus, p) ...
        - canonicalDynamics(minus, p)) / (2 * step);
end
assert(max(abs(analyticJacobian - finiteJacobian), [], 'all') < 1e-6);

% Verify the Hamiltonian gradient against the same central-difference scale.
analyticGradient = hamiltonianGradient(canonical, p);
finiteGradient = zeros(14, 1);
for column = 1:14
    step = 1e-6 * max(1, abs(canonical(column)));
    plus = canonical;
    minus = canonical;
    plus(column) = plus(column) + step;
    minus(column) = minus(column) - step;
    finiteGradient(column) = (hamiltonian(plus, p) ...
        - hamiltonian(minus, p)) / (2 * step);
end
assert(max(abs(analyticGradient - finiteGradient)) < 1e-6);

% Verify that identical inertial histories produce exactly zero NTC errors.
reference = [53.15, 0, 0, 0, 2.8e-3, 1.7e-4; ...
    53.14, 0.1, 0.01, -5e-6, 2.8e-3, 1.7e-4];
[positionNtc, velocityNtc] = ntcDifferences(reference, reference);
assert(all(positionNtc == 0, 'all'));
assert(all(velocityNtc == 0, 'all'));

fprintf('Derivative and NTC verification passed.\n');
