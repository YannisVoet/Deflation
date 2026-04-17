function[T]=def_param(precond_tech)

% DEF_PARAM: Set default parameters for preconditioning techniques.
%
% T = DEF_PARAM(prec) initializes the default parameters for a custom list
% of preconditioners. Available preconditioners are:
% 'No preconditioning'  - identity.
% 'Jacobi'              - Jacobi preconditioner (or diagonal scaling)
% 'Cholesky'            - Incomplete Cholesky preconditioner
% 'SIPIC'               - SIPIC preconditioner (see [1])
% 'Schwarz'             - Additive Schwarz preconditioner (see [2])
% 'Deflation'           - Deflation-based preconditioner (see [3])
%
% References:
% [1] F. de Prenter, C. V. Verhoosel, G. J. van Zwieten, and E. H. van Brummelen. 
% Condition number analysis and preconditioning of the finite cell method. 
% Computer Methods in Applied Mechanics and Engineering, 316:297–327, 2017.
% [2] F. de Prenter, C. Verhoosel, and E. Van Brummelen. Preconditioning 
% immersed isogeometric finite element methods with application to flow problems. 
% Computer Methods in Applied Mechanics and Engineering, 348:604–631, 2019.
% [3] Y. Voet, M. Möller, P. Antolin and K. Vuik. Deflation-based preconditioning 
% for immersed finite element methods and immersogeometric analysis, 2026.

% Preconditioning techniques
def_tech={'No preconditioning', 'Jacobi', 'Cholesky', 'SIPIC', 'Schwarz', 'Deflation'};
colors={[0 0 0], [0 0.2 0.8], [0.7 0 1], [0.4 0.7 0.1], [0.8 0.8 0.2], [0.9 0.1 0.1]};
markers={'o', 'x', '+', '*', '^', 's'};

% Deflation default parameters
param_defl.gamma=1; % Volume threshold for identifying elements with small volume ratio
param_defl.reduct=false; % Activate deflation rank-reduction algorithm

% Schwarz default parameters
param_sch.block_sel='cut_elem'; % One index block for each cut element
param_sch.inv='exact'; % Approximate inverse for improved stability
param_sch.tol=1e-14; % Truncation tolerance

% Incomplete Cholesky default parameters
param_chol.type='nofill';

params={[],[],param_chol,[],param_sch,param_defl};
n_prec=length(def_tech);

for k=1:n_prec
    T(k).name=def_tech{k};
    T(k).color=colors{k};
    T(k).marker=markers{k};
    T(k).addendum='';
    T(k).param=params{k};
end

[~,ids]=ismember(precond_tech, def_tech);
T=T(ids);
end