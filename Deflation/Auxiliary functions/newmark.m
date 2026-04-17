function[u, varargout]=newmark(M, C, K, f, u0, varargin)

% NEWMARK: Numerical approximation of the solution to the system of
% differential equations M*u'' + C*u' + K*u = f using the Newmark method.
%
% u = NEWMARK(M, C, K, f, u0) computes the numerical solution on the 
% interval I = [0 1] using the implicit unconditionally stable version of 
% the Newmark method (with parameters [beta, gamma] = [1/4 1/2]) and 
% N = 200 subdivisions in time. The matrices M, C and K are n x n matrices, 
% f is a function handle such that f(t) returns a vector of size n. The 
% initial conditions are stored in the n x 1 x 2 tensor u0. The solution u 
% is a n x (N + 1) x 3 tensor storing the solution and its first and second
% time derivatives along the first, second and third pages, respectively.
%
% u = NEWMARK(M, C, K, f, u0, Options) specifies an Options structure 
% (see the list of available options below).
%
% u = NEWMARK(M, C, K, f, u0, name, value) specifies optional name/value 
% pairs. Available options are:
% 'I'       - computes the solution on an interval I. Default: I = [0 1].
% 'N'       - specifies the number of subdisions. Default: N = 200.
% 'param'   - specifies the [beta, gamma] parameters of the method.
%             Default: [1/4 1/2] (implicit unconditionally stable version).
%
% [u, Output] = NEWMARK(M, C, K, f, u0, ...) also returns a structure with 
% diagnostic information, including the computing time for solving linear 
% systems.

%% Set algorithm parameters

n=size(M,1);

Param = inputParser;
Param.KeepUnmatched = true;
Param.addRequired('M', @(x) all(size(x)==[n n]));
Param.addRequired('C', @(x) all(size(x)==[n n]));
Param.addRequired('K', @(x) all(size(x)==[n n]));
Param.addRequired('f', @(x) isa(x, 'function_handle'));
Param.addRequired('u0', @(x) all(size(x)==[n 1 2]));
Param.addParameter('I', [0 1], @(x) isa(x, 'double') && length(x)==2);
Param.addParameter('N', 200, @(x) isscalar(x) & x > 0);
Param.addParameter('param', [1/4 1/2], @(x) isa(x, 'double') && length(x)==2);
Param.parse(M, C, K, f, u0, varargin{:});

%% Retrieve parameters
I=Param.Results.I;
N=Param.Results.N;
param=Param.Results.param;

%% Pre-processing
Dt=(I(2)-I(1))/N;
T=linspace(I(1), I(2), N+1);

%% Newmark method

% (beta,gamma) parameters
beta=param(1);
gamma=param(2);

A=M+gamma*Dt*C+beta*Dt^2*K;

% Sparse direct solver
dA=decomposition(A);

% Initialization
u=zeros(n,N+1,3);
u(:,1,1:2)=u0;
u(:,1,3)=M\(f(T(1))-K*u(:,1,1)-C*u(:,1,2));

% Measure the time for solving linear systems
time=0;

for n=1:N

    % Computation of the solution
    ut=u(:,n,1)+Dt*u(:,n,2)+(0.5-beta)*Dt^2*u(:,n,3);
    vt=u(:,n,2)+(1-gamma)*Dt*u(:,n,3);
    rhs=f(T(n+1))-C*vt-K*ut;
    tic
    u(:,n+1,3)=dA\rhs;
    time=time+toc;
    u(:,n+1,2)=vt+gamma*Dt*u(:,n+1,3);
    u(:,n+1,1)=ut+beta*Dt^2*u(:,n+1,3);

end

% Create output structure
Output.time=time;
varargout{1}=Output;
end