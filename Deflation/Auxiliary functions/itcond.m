function[kappa, min_mod, max_mod] = itcond(A, k, sigma, B)

% ITCOND: Computes the iterative condition number of A, defined as the
% ratio of largest to smallest moduli of eigenvalues.
%
% kappa = ITCOND(A) returns the iterative condition number.
% 
% kappa = ITCOND(A, s) ignores the s eigenvalues of smallest moduli.
% Default: 0.
%
% kappa = ITCOND(A, s, sigma) computes the smallest eigenvalues closest to a
% shift sigma. Default: 0.
%
% kappa = ITCOND(A, s, sigma, B) computes the iterative condition number of 
% the matrix pair (A,B), where B is assumed positive definite.
%
% [kappa, min_mod, max_mod] = ITCOND(A, s, sigma, B) also returns the 
% eigenvalues of smallest and largest modulus.

arguments
    A
    k=0
    sigma=0;
    B=[];
end

if isempty(k)
    k=0;
end

if isempty(B)
    max_mod=abs(eigs(A, 1, 'largestabs'));
    min_mod=max(abs(eigs(A, k+1, sigma, 'SubspaceDimension', min([100, size(A,1)]))));
else
    max_mod=abs(eigs(A, B, 1, 'largestabs'));
    min_mod=max(abs(eigs(A, B, k+1, sigma, 'SubspaceDimension', min([100, size(A,1)]))));
end

kappa = max_mod/min_mod;

if ~isreal(kappa)
    disp('Eigenvalues failed to converge for condition number computation.')
end
end