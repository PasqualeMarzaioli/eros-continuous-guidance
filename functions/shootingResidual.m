%SHOOTINGRESIDUAL  Terminal PMP residual and analytical shooting Jacobian.
%
%   Propagates the 14-D canonical system (and costate sensitivities when
%   gradients are requested) and returns boundaryResidual at tf.
%
%   Author: Pasquale Marzaioli

function [residual, jacobian] = shootingResidual(decision, xInitial, ...
        target, p, odeOptions)
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
