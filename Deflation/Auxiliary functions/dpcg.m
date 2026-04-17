function[x, varargout] = dpcg(A, b, varargin)

% DPCG: Deflated Preconditioned Conjugate Gradient algorithm for solving
% linear systems of equations. Based on Algorithm 6.18 of [1] and
% Algorithm 1 of [2].
%
% x = DPCG(A, b) iteratively solves the linear system A*x = b where A is a
% n x n symmetric positive definite matrix and b is a column vector of
% length n. The matrix A is either provided explicitly or as a function
% handle.
%
% x = DPCG(A, b, tol) specifies the tolerance on the absolute residual. For a
% stopping criterion based on the relative residual, replace tol with
% tol*norm(b). Default: 1e-8.
%
% x = DPCG(A, b, tol, maxiter) also specifies the maximum number of
% iterations. Default: min(n,10).
%
% x = DPCG(A, b, tol, maxiter, P) also specifies a projection matrix for
% deflation. The projection matrix P is either provided explicitly or as a
% function handle. Note: PA must be symmetric. Default: identity.
%
% x = DPCG(A, b, tol, maxiter, P, M) also specifies a preconditioning matrix.
% The preconditioning matrix M is either provided explicitly or as a
% function handle. Default: identity.
%
% x = DPCG(A, b, tol, maxiter, P, M, x0) also specifies the initial starting
% vector. Default: all zero vector.
%
% x = DPCG(A, b, tol, maxiter, P, M, x0, S) also specifies the "uniqueness 
% operator" for recovering the solution of Ax = b. See Algorithm 1 in [3].
%
% x = DPCG(A, b, tol, maxiter, P, M, x0, S, name, value) also specifies
% optional name/value pair arguments. Available parameters are:
% 'xt'      - "true" solution of the linear system (used for convergence analyses).
% 'store'   - Boolean indicator for storing intermediate solutions. If
%             'store' is true, x is a matrix containing the approximate
%             solutions at each step along its columns. Default: false.
%
% [x, niter] = DPCG(A, b, ...) also returns the number of iterations.
%
% [x, niter, res] = DPCG(A, b, ...) returns the vector of absolute residuals 
% for each iteration.
%
% [x, niter, res, res_prec] = DPCG(A, b, ...) returns the vector of absolute 
% preconditioned residuals for each iteration.
%
% [x, niter, res, res_prec, err] = DPCG(A, b, ...) returns the vector of errors 
% for each iteration (if the true solution was specified).
%
% Reference:
% [1] Y. Saad. Iterative methods for sparse linear systems. SIAM, 2003.
% [2] F. Vermolen, K. Vuik, and G. Segal. Deflation in preconditioned
% conjugate gradient methods for finite element problems. In Conjugate
% Gradient Algorithms and Finite Element Methods, Springer, 2004.
% [3] J. M. Tang, R. Nabben, C. Vuik, and Y. A. Erlangga. Comparison of 
% two-level preconditioners derived from deflation, domain decomposition 
% and multigrid methods. Journal of scientific computing, 39(3):340–370, 2009.
%
% See also DEFLATION, CG, GMRESK, BICGSTB.

%% Set algorithm parameters

% Set default parameters
Default{1}=1e-8;
Default{2}=min(length(b),10);
Default{3}=@(x) x;
Default{4}=@(x) x;
Default{5}=zeros(size(b));
Default{6}=@(x) x;
Default{7}=zeros(size(b));
Default{8}=false;

% Replace empty inputs with default parameters
def=cell2mat(cellfun(@isempty, varargin, 'UniformOutput', false));
[varargin{def}]=Default{def};

Param = inputParser;
Param.KeepUnmatched = true;
Param.addRequired('A', @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addRequired('b');
Param.addOptional('tol',        Default{1}, @(x) isscalar(x) && x > 0);
Param.addOptional('maxiter',    Default{2}, @(x) isscalar(x) && x > 0);
Param.addOptional('P',          Default{3}, @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addOptional('M',          Default{4}, @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addOptional('x0',         Default{5}, @(x) all(size(x)==size(b)));
Param.addOptional('S',          Default{6}, @(x) isa(x, 'function_handle'));
Param.addParameter('xt',        Default{7}, @(x) all(size(x)==size(b)));
Param.addParameter('store',     Default{8}, @(x) islogical(x));
Param.parse(A, b, varargin{:});

%% Retrieve parameters
tol=Param.Results.tol;
maxiter=Param.Results.maxiter;
P=Param.Results.P;
M=Param.Results.M;
x=Param.Results.x0;
xt=Param.Results.xt;
S=Param.Results.S;
store=Param.Results.store;

if isa(A, 'double')
    A=@(x) A*x;
end

if isa(P, 'double')
    P=@(x) P*x;
end

if isa(M, 'double')
    dM=decomposition(M);
    M=@(x) dM\x;
end

if store % store all intermediate solutions
    X=cell(1,maxiter+1);
    X{1}=x;
end

%% DPCG method

r=P(b-A(x));
z=M(r);
p=z;
q=P(A(p));
xi=dot(p,q);
eta=dot(r,z);

res=zeros(maxiter+1,1);
res_prec=zeros(maxiter+1,1);
err=zeros(maxiter+1,1);

res(1)=norm(r);
res_prec(1)=sqrt(eta);
err(1)=sqrt((xt-S(x))'*A(xt-S(x)));

for k=1:maxiter

    alpha=eta/xi;
    x=x+alpha*p;
    r=r-alpha*q;
    z=M(r);
    eta=dot(r,z);
    beta=eta/(alpha*xi);
    p=z+beta*p;
    q=P(A(p));
    xi=dot(p,q);
    res(k+1)=norm(r);
    res_prec(k+1)=sqrt(eta);
    err(k+1)=sqrt((xt-S(x))'*A(xt-S(x)));

    if store
        X{k+1}=x;
    end

    if res_prec(k+1)<tol
        break
    end

end

if store
    x=cell2mat(X(1:k+1));
end

x=S(x);
varargout{1}=k;
varargout{2}=res(1:k+1);
varargout{3}=res_prec(1:k+1);
varargout{4}=err(1:k+1);

end