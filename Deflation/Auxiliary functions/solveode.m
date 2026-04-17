function [Data, Solution, varargout] = solveode(~, Data, M, C, K, varargin)

% SOLVEODE: Numerical approximation of the solution to the system of
% differential equations M*u'' + C*u' + K*u = f.
%
% [Data, Solution] = SOLVEODE(Mesh, Data, M, C, K) computes the numerical
% solution as well as its first and second derivatives on the interval
% I = [0 1] using N = 200 subdivisions in time. The parameters of the
% method are stored in the Data structure and the numerical solution is
% stored in the solution structure. The matrices M, C and K are the global
% mass, damping and stiffness matrices, respectively. The right-hand side
% function and Neumann boundary conditions may be specified as a function
% handle or as a cell array {{f1x, f2x, ..., fkx}, ft} containing the
% spatial and temporal components of separable functions such that
% f = ∑ fxi(x)*fti(t). Time dependency of data is specified using the
% variable 't'.
%
% [Data, Solution] = SOLVEODE(Mesh, Data, M, C, K, Options) specifies an
% Options structure (see the list of available options below).
%
% [Data, Solution] = SOLVEODE(Mesh, Data, M, C, K, name, value) specifies
% options as name/value pairs. Available options are:
% 'I'       - computes the solution on an interval I. Default: I = [0 1].
% 'N'       - specifies the number of subdisions. Default: N = 200.
% 'V'       - specifies a reduced basis and implicitly solves
%             Mt*x'' + Ct*x' + Kt*x = ft where u = V*x, ft = V'*f and
%             Mt = V'*M*V, Ct = V'*C*V, Kt = V'*K*V.
% 'S'       - specifies a change of basis and implicitly solves
%             Mt*x'' + Ct*x' + Kt*x = ft where u = S*x, ft = S'*f and
%             Mt = S'*M*S, Ct = S'*C*S, Kt = S'*K*S.
%             Note: the change of basis is performed on the global system,
%             contrary to the reduced basis approach.
% 'Sinv'    - Inverse of the change of basis.
% 'M11'     - mass matrix provided to the time-stepping scheme.
%             Example: M11 is a lumped mass matrix.
% 'C11'     - damping matrix provided to the time-stepping scheme.
% 'K11'     - stiffness matrix provided to the time-stepping scheme.
% 'solver'  - specifies the numerical time integration scheme. Available
%             schemes are 'newmark' (default) and 'erk'.
% 'lumping' - specifies a lumping scheme to apply to the global mass matrix
%             as a function handle. Default: id (i.e. no lumping).
%             Example: @(x) lump(x, 1);
% 'param'   - specifies optional parameters for the time integration
%             scheme (see the options supported for each solver).
%
% [Data, Solution, Output] = SOLVEODE(Mesh, Data, M, C, K, ...) also returns
% a structure with diagnostic information, including the computing time
% for solving linear systems.

% SOLVEODE(Coeff, Data, Options)
% Coeff: struct array with fields:
%   Ai  - ith coefficient matrix for the ODE
%         An*u^(n) + ... + A1*u^(1) + A0*u^(0) = f.
%         Solvers currently supported (2nd order ODEs)
% Data: struct array with fields:

%% Set algorithm parameters

n=size(M,1);

Param = inputParser;
Param.KeepUnmatched = true;
Param.addRequired('Data', @isstruct);
Param.addRequired('M', @(x) all(size(x)==[n n]));
Param.addRequired('C', @(x) all(size(x)==[n n]));
Param.addRequired('K', @(x) all(size(x)==[n n]));
Param.addParameter('I', [0 1], @(x) isa(x, 'double') && length(x)==2);
Param.addParameter('N', 200, @(x) isscalar(x) & x > 0);
Param.addParameter('V', [], @ismatrix);
Param.addParameter('S', speye(n), @ismatrix);
Param.addParameter('Sinv', speye(n), @ismatrix);
Param.addParameter('M11', [], @ismatrix);
Param.addParameter('C11', [], @ismatrix);
Param.addParameter('K11', [], @ismatrix);
Param.addParameter('solver', 'newmark', @(x) ismember(x,{'newmark','erk'}));
Param.addParameter('lumping', @(x) x, @(x) isa(x, 'function_handle'));
Param.parse(Data, M, C, K, varargin{:});

%% Retrieve parameters
I=Param.Results.I;
N=Param.Results.N;
V=Param.Results.V;
S=Param.Results.S;
Sinv=Param.Results.Sinv;
solver=Param.Results.solver;
lumping=Param.Results.lumping;

%% Pre-processing
delta_t=(I(2)-I(1))/N;
T=linspace(I(1), I(2), N+1);

% Saving data
Data.Discretization.N=N;
Data.Discretization.I=I;
Data.Discretization.delta_t=delta_t;
Data.Discretization.time_vec=T;

% Check time dependency of data
% 'rhs', 'nmnn' and 'drchlt' are boolean indicators for time-dependent
% data. The algorithm may be significantly faster if data is
% time-independent.
if iscell(Data.f)
    rhs=dependency(Data.f{2},'t');
else
    rhs=dependency(Data.f,'t');
end

if iscell(Data.h)
    nmnn=dependency(Data.h{2},'t');
else
    nmnn=dependency(Data.h,'t');
end

if iscell(Data.g)
    drchlt=dependency(Data.g{2},'t');
else
    drchlt=dependency(Data.g,'t');
end


% Time dependency
Data.rhs=rhs;
Data.drchlt=drchlt;
Data.nmnn=nmnn;

% Reduced basis
mor=~isempty(V);

% Retrieve data parameters
% Interior degrees of freedom
Nf=Data.Nf;
% Dirichlet degrees of freedom
Nd=Data.Nd;

%% Compute initial conditions
[u]=initcond(Data);

%% Compute Dirichlet boundary conditions

[ud,Data]=dirichlet(Data, T(1));

% Evaluation
if drchlt
    for n=1:N+1
        [ud]=dirichlet(Data, T(n));
        u(Nd,n,:)=ud(Nd,:,:);
    end
else
    u(Nd,:,:)=repmat(ud(Nd,:,:), [1 N+1]);
end

%% Time integration

% Change of basis
u=sppagemtimes(Sinv,u);

M=S'*M*S;
C=S'*C*S;
K=S'*K*S;

M=lumping(M);

Data.S=S;
Data.Sinv=Sinv;

% Partitioning of M, C and K in block matrices
M11=M(Nf,Nf);
C11=C(Nf,Nf);
K11=K(Nf,Nf);

Data.M12=M(Nf,Nd);
Data.C12=C(Nf,Nd);
Data.K12=K(Nf,Nd);

[~,Data]=eval_rhs(Data, T(1), 'init', varargin{:});

if mor
    m=size(V,2);
    x0=zeros(m,1,2);

    % Least squares minimizer in the norm induced by M11
    x0(:,1)=V'*(M11*u(Nf,1,1));
    x0(:,2)=V'*(M11*u(Nf,1,2));

    M11=V'*M11*V;
    K11=V'*K11*V;
    C11=V'*C11*V;

    x0(:,1)=M11\x0(:,1);
    x0(:,2)=M11\x0(:,2);

    [x,Output]=feval(solver, coeff(M11,'M11',Param), coeff(C11,'C11',Param), coeff(K11,'K11',Param), @(t) V'*eval_rhs(Data, t, 'eval', varargin{:}), x0, varargin{:});

    u(Nf,:,:)=sppagemtimes(V,x);
else
    [u(Nf,:,:),Output]=feval(solver, coeff(M11,'M11',Param), coeff(C11,'C11',Param), coeff(K11,'K11',Param), @(t) eval_rhs(Data, t, 'eval', varargin{:}), u(Nf,1,1:2), varargin{:});
    u=sppagemtimes(S,u);
end

% Create solution structure
Solution.U=u(:,:,1);
Solution.d1_U=u(:,:,2);
Solution.d2_U=u(:,:,3);

% Return output structure
varargout{1}=Output;

    function[M]=coeff(M, name, Param)
        if ~isempty(Param.Results.(name))
            M=Param.Results.(name);
        end
    end
end