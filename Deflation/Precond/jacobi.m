function[Ap, D] = jacobi(A)

% JACOBI: Computes the diagonally preconditioned matrix D*A*D, where
% D = 1./sqrt(diag(A)). The matrix A must be real positive definite (but
% not necessarily symmetric). If A is the Gram matrix for a bilinear form
% a : Vh x Vh -> R with respect to a basis Phi = {phi_1,...,phi_n},
% diagonal preconditioning amounts to forming the Gram matrix in the
% rescaled basis Phih = {phih_1,...,phih_n}, where phih_j =
% phi_j/||phi_j||_a and ||u||_a = sqrt(a(u,u)).
%
% Ap = JACOBI(A) returns the preconditioned matrix Ap = D*A*D.
%
% [Ap, D] = JACOBI(A) also returns the diagonal matrix D.

% Check if A is square
[m, n] = size(A);
if m ~= n
    error('Matrix A must be square.');
end

% Extract diagonal elements
diag_A = diag(A);

% Check if rescaling is possible
if any(diag_A <= 0)
    error('The matrix is not positive definite. Rescaling is impossible.');
end

% Compute diagonal scaling matrix
D = spdiags(1 ./ sqrt(diag_A), 0, n, n);

% Apply symmetric diagonal scaling
Ap = D * A * D;

end