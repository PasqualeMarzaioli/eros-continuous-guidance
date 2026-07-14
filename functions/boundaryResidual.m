%BOUNDARYRESIDUAL  Terminal position/velocity, free-mass, and free-time conditions.
%
%   Eight conditions: r(tf)=r_f, v(tf)=v_f, lambda_m(tf)=0 (free final mass),
%   and H(tf)=0 (free final time transversality for the autonomous problem).
%
%   Author: Pasquale Marzaioli

function [residual, stateJacobian] = boundaryResidual(canonical, target, p)
residual = [canonical(1:6) - target; canonical(14); ...
    hamiltonian(canonical, p)];

if nargout > 1
    stateJacobian = zeros(8, 14);
    stateJacobian(1:6, 1:6) = eye(6);
    stateJacobian(7, 14) = 1;
    stateJacobian(8, :) = hamiltonianGradient(canonical, p).';
end
end
