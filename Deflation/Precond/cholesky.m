function[Ap, L, Linv] = cholesky(A, param)

% CHOLESKY: Computes the incomplete Cholesky preconditioned matrix
% L^{-1}*A*L^{-T}, where L is the incomplete Cholesky factor of A. The
% matrix A must be real positive definite (but not necessarily symmetric).
%
% Ap = CHOLESKY(A) returns the preconditioned matrix Ac = L^{-1}*A*L^{-T}.
%
% [Ap, L] = CHOLESKY(A) also returns the incomplete Cholesky factor of A.
%
% [Ap, L, Linv] = CHOLESKY(A) also returns the inverse of the incomplete 
% Cholesky factor of A.
%
% Ap = CHOLESKY(A, param) specifies a structure of parameters for the
% incomplete Cholesky factorization. Those are the same parameters as
% MATLAB's built-in ichol function.

% Check if A is square
[m, n] = size(A);
if m ~= n
    error('Matrix A must be square.');
end

% Incomplete Cholesky factorization
if nargin > 1
    fields=fieldnames(param);
    rmfields=setdiff(fields, {'type', 'droptol', 'michol', 'diagcomp', 'shape'});
    param=rmfield(param, rmfields);
    L = ichol(A, param);
else
    L = ichol(A);
end

% Apply symmetric diagonal scaling
Ap = L\(A/L');

if nargout > 2
    % Inverse Choleky factor
    Linv = inv(L);
end
end