function[x, varargout] = cg(A, b, varargin)

% CG: Conjugate Gradient algorithm for solving linear systems of equations.
% Based on Algorithm 6.18 of [1].
%
% x = CG(A, b) iteratively solves the linear system A*x = b where A is a
% n x n symmetric positive definite matrix and b is a column vector of
% length n. The matrix A is either provided explicitly or as a function
% handle.
%
% x = CG(A, b, tol) specifies the tolerance on the absolute residual. For a
% stopping criterion based on the relative residual, replace tol with
% tol*norm(b). Default: 1e-8.
%
% x = CG(A, b, tol, maxiter) also specifies the maximum number of
% iterations. Default: min(n,10).
%
% x = CG(A, b, tol, maxiter, M) also specifies a preconditioning matrix.
% The preconditioning matrix M is either provided explicitly or as a
% function handle. Default: identity.
%
% x = CG(A, b, tol, maxiter, M, x0) also specifies the initial starting
% vector. Default: all zero vector.
%
% x = CG(A, b, tol, maxiter, M, x0, name, value) also specifies
% optional name/value pair arguments. Available parameters are:
% 'xt'      - "true" solution of the linear system (used for convergence analyses).
% 'store'   - Boolean indicator for storing intermediate solutions. If
%             'store' is true, x is a matrix containing the approximate
%             solutions at each step along its columns. Default: false.
%
% [x, niter] = CG(A, b, ...) also returns the number of iterations.
%
% [x, niter, res, res_prec, err] = CG(A, b, ...) also returns the vector of
% absolute residuals res, absolute preconditioned residuals res_prec and 
% errors err (if xt is supplied) for each iteration.
%
% Reference:
% [1] Y. Saad. Iterative methods for sparse linear systems. SIAM, 2003.
%
% See also GMRESK, BICGSTB.

%% Set algorithm parameters

% Set default parameters
Default{1}=1e-8;
Default{2}=min(length(b),10);
Default{3}=@(x) x;
Default{4}=zeros(size(b));
Default{5}=zeros(size(b));
Default{6}=false;

% Replace empty inputs with default parameters
def=cell2mat(cellfun(@isempty, varargin, 'UniformOutput', false));
[varargin{def}]=Default{def};

Param = inputParser;
Param.KeepUnmatched = true;
Param.addRequired('A', @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addRequired('b');
Param.addOptional('tol',        Default{1}, @(x) isscalar(x) && x > 0);
Param.addOptional('maxiter',    Default{2}, @(x) isscalar(x) && x > 0);
Param.addOptional('M',          Default{3}, @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addOptional('x0',         Default{4}, @(x) all(size(x)==size(b)));
Param.addParameter('xt',        Default{5}, @(x) all(size(x)==size(b)));
Param.addParameter('store',     Default{6}, @(x) islogical(x));
Param.parse(A, b, varargin{:});

%% Retrieve parameters
tol=Param.Results.tol;
maxiter=Param.Results.maxiter;
M=Param.Results.M;
x=Param.Results.x0;
xt=Param.Results.xt;
store=Param.Results.store;

if isa(A, 'double')
    A=@(x) A*x;
end

if isa(M, 'double')
    dM=decomposition(M);
    M=@(x) dM\x;
end

if store % store all intermediate solutions
    X=cell(1,maxiter+1);
    X{1}=x;
end

%% CG method

r=b-A(x);
z=M(r);
p=z;
q=A(p);
xi=dot(p,q);
eta=dot(r,z);

res=zeros(maxiter+1,1);
res_prec=zeros(maxiter+1,1);
err=zeros(maxiter+1,1);

res(1)=norm(r);
res_prec(1)=sqrt(eta);
err(1)=sqrt((xt-x)'*A(xt-x));


for k=1:maxiter

    alpha=eta/xi;
    x=x+alpha*p;
    r=r-alpha*q;
    z=M(r);
    eta=dot(r,z);
    beta=eta/(alpha*xi);
    p=z+beta*p;
    q=A(p);
    xi=dot(p,q);
    res(k+1)=norm(r);
    res_prec(k+1)=sqrt(eta);
    err(k+1)=sqrt((xt-x)'*A(xt-x));

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

varargout{1}=k;
varargout{2}=res(1:k+1);
varargout{3}=res_prec(1:k+1);
varargout{4}=err(1:k+1);

end