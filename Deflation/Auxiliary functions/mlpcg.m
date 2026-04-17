function[x, varargout] = mlpcg(A, b, varargin)

% MLPCG: Multi-level Preconditioned Conjugate Gradient algorithm for solving
% linear systems of equations. Based on Algorithm 1 of [1].
%
% x = MLPCG(A, b) iteratively solves the linear system A*x = b where A is a
% n x n symmetric positive definite matrix and b is a column vector of
% length n. The matrix A is either provided explicitly or as a function
% handle.
%
% x = MLPCG(A, b, tol) specifies the tolerance on the absolute residual. For a
% stopping criterion based on the relative residual, replace tol with
% tol*norm(b). Default: 1e-8.
%
% x = MLPCG(A, b, tol, maxiter) also specifies the maximum number of
% iterations. Default: min(n,10).
%
% x = MLPCG(A, b, tol, maxiter, M1, M2, M3) also specifies the operators M1,
% M2 and M3 either provided explicitly or as function handles. Default: identity.
%
% x = MLPCG(A, b, tol, maxiter, M1, M2, M3, x0) also specifies the initial 
% starting vector. Default: all zero vector.
%
% x = MLPCG(A, b, tol, maxiter, M1, M2, M3, x0, name, value) also specifies
% optional name/value pair arguments. Available parameters are:
% 'Vs'      - starting vector.
% 'Ve'      - ending vector ("uniqueness step"). See [1].
% 'xt'      - "true" solution of the linear system (used for convergence analyses).
% 'store'   - Boolean indicator for storing intermediate solutions. If
%             'store' is true, x is a matrix containing the approximate
%             solutions at each step along its columns. Default: false.
%
% [x, niter] = MLPCG(A, b, ...) also returns the number of iterations.
%
% [x, niter, res, res_prec, err] = MLPCG(A, b, ...) also returns the vector of
% absolute residuals res, absolute preconditioned residuals res_prec and 
% errors err (if xt is supplied) for each iteration.
%
% Reference:
% [1] J. M. Tang, R. Nabben, C. Vuik, and Y. A. Erlangga. Comparison of 
% two-level preconditioners derived from deflation, domain decomposition 
% and multigrid methods. Journal of scientific computing, 39(3):340–370, 2009.
%
% See also CG, GMRESK, BICGSTB.

%% Set algorithm parameters

% Set default parameters
Default{1}=1e-8;
Default{2}=min(length(b),10);
Default{3}=@(x) x;
Default{4}=@(x) x;
Default{5}=@(x) x;
Default{6}=zeros(size(b));
Default{7}=zeros(size(b));
Default{8}=@(x) x;
Default{9}=@(x) x;
Default{10}=false;

% Replace empty inputs with default parameters
def=cell2mat(cellfun(@isempty, varargin, 'UniformOutput', false));
[varargin{def}]=Default{def};

Param = inputParser;
Param.KeepUnmatched = true;
Param.addRequired('A', @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addRequired('b');
Param.addOptional('tol',        Default{1}, @(x) isscalar(x) && x > 0);
Param.addOptional('maxiter',    Default{2}, @(x) isscalar(x) && x > 0);
Param.addOptional('M1',         Default{3}, @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addOptional('M2',         Default{4}, @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addOptional('M3',         Default{5}, @(x) isa(x, 'function_handle') || isa(x, 'double'));
Param.addOptional('x0',         Default{6}, @(x) all(size(x)==size(b)));
Param.addParameter('xt',        Default{7}, @(x) all(size(x)==size(b)));
Param.addParameter('Vs',        Default{8}, @(x) isa(x, 'function_handle'));
Param.addParameter('Ve',        Default{9}, @(x) isa(x, 'function_handle'));
Param.addParameter('store',     Default{10}, @(x) islogical(x));
Param.parse(A, b, varargin{:});

%% Retrieve parameters
tol=Param.Results.tol;
maxiter=Param.Results.maxiter;
M1=Param.Results.M1;
M2=Param.Results.M2;
M3=Param.Results.M3;
x=Param.Results.x0;
xt=Param.Results.xt;
Vs=Param.Results.Vs;
Ve=Param.Results.Ve;
store=Param.Results.store;

if isa(A, 'double')
    A=@(x) A*x;
end

if isa(M1, 'double')
    M1=@(x) M1*x;
end

if isa(M2, 'double')
    M2=@(x) M2*x;
end

if isa(M3, 'double')
    M3=@(x) M3*x;
end

if store % store all intermediate solutions
    X=cell(1,maxiter+1);
    X{1}=x;
end

%% PCG method

x=Vs(x);
% r=b-A(x);
r=M3(b-A(x));
y=M1(r);
p=M2(y);
w=M3(A(p));
xi=dot(p,w);
eta=dot(r,y);

res=zeros(maxiter+1,1);
res_prec=zeros(maxiter+1,1);
err=zeros(maxiter+1,1);

res(1)=norm(r);
res_prec(1)=sqrt(eta);
err(1)=sqrt((xt-Ve(x))'*A(xt-Ve(x)));

for k=1:maxiter

    alpha=eta/xi;
    x=x+alpha*p;
    r=r-alpha*w;
    y=M1(r);
    eta=dot(r,y);
    beta=eta/(alpha*xi);
    p=M2(y)+beta*p;
    w=M3(A(p));
    xi=dot(p,w);
    res(k+1)=norm(r);
    res_prec(k+1)=sqrt(eta);
    err(k+1)=sqrt((xt-Ve(x))'*A(xt-Ve(x)));

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

x=Ve(x);
varargout{1}=k;
varargout{2}=res(1:k+1);
varargout{3}=res_prec(1:k+1);
varargout{4}=err(1:k+1);

end