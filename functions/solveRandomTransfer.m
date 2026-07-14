%SOLVERANDOMTRANSFER  Random costate screening then PMP single-shooting refinement.
%
%   Draws batches of initial costates (wide and order-one scales, planar
%   lambda_z = lambda_vz = 0), ranks by cheap residual, and refines the
%   best candidates from every configured batch; keeps the lowest-cost
%   converged extremal found by this finite multistart search.
%
%   Author: Pasquale Marzaioli

function [solution, time, trajectory] = solveRandomTransfer(xInitial, target, ...
        p, odeOptions, screeningOptions, candidateCount, maximumBatches, seed)
rng(seed, 'twister');
initialPeriod = 2 * pi * norm(xInitial(1:3))^(3 / 2);
bestSolution = [];
bestCost = inf;
bestBatch = NaN;

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
                target, p, odeOptions, false);
            [trialTime, trialTrajectory] = propagateCanonical(...
                trial, xInitial, p, odeOptions, 1200);
            trialResidual = shootingResidual(trial, xInitial, target, ...
                p, odeOptions);
            trialCost = trapz(trialTime, dustDensity(...
                vecnorm(trialTrajectory(:, 1:3), 2, 2), p));

            if norm(trialResidual) < 1e-8 && trialCost < bestCost
                bestCost = trialCost;
                bestSolution = trial;
                bestBatch = batch;
            end
        catch
            % A failed random candidate is discarded; the next candidate is independent.
        end
    end

end

if isempty(bestSolution)
    error('No converged PMP extremal was found after %d random batches.', ...
        maximumBatches);
end

fprintf('Selected batch %d anchor with nondimensional exposure %.10e.\n', ...
    bestBatch, bestCost);

solution = bestSolution;
[time, trajectory] = propagateCanonical(...
    solution, xInitial, p, odeOptions, 1800);
end
