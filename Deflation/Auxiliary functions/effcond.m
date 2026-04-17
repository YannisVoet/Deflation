function[kappa] = effcond(A,s)

% EFFCOND: Computes the (effective) condition number of the m x n matrix A.
%
% kappa = EFFCOND(A) computes the condition number sigma_1(A)/sigma_l(A),
% where l = min(m,n).
%
% kappa = EFFCOND(A,s) computes the effective condition number of A, defined as
% sigma_1(A)/sigma_r(A), where r is the rank of the m x n matrix A and s is 
% the number of zero singular values (s = min(m,n) - r).

arguments
    A
    s=0;
end

sigma1=svds(A, 1, 'largest');
sigmar=max(svds(A, s+1, 'smallest'));
kappa=sigma1/sigmar;

if ~isreal(kappa)
    disp('Singular values failed to converge for condition number computation.')
end
end