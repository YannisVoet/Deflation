function[nr]=normA(A,x)

% NORMA: Computes the A-norm of x: ||x||_A = sqrt(x'*A*x) for a Hermitian
% positive definite matrix A.
%
% z = NORMA(A,x) returns the A-norm of x. A may either be provided
% explicitly or as a function handle.

if isa(A, 'double')
    A=@(x) A*x;
end

nr=sqrt(x'*A(x));
end