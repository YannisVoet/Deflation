function[Ap,P1,P2,Pfun,Ptfun,r,iss,isl,Q]=deflation(sp_trimmed, msh_trimmed, A, int_dofs, varargin)

% DEFLATION: Computes the deflation-based preconditioner.
%
% Ap = DEFLATION(sp_trimmed, msh_trimmed, A, int_dofs) returns the
% preconditioned matrix Ap = P*D*A*D for the composition of diagonal scaling 
% and deflation preconditioning, where D = 1./sqrt(diag(A)) is the Jacobi 
% preconditioner and P is the deflation preconditioner.
%
% Ap = DEFLATION(sp_trimmed, msh_trimmed, A, int_dofs, name, value) specifies 
% optional name/value pair arguments. Available parameters are:
% 'gamma'       - only deflates the elements whose volume fraction is smaller 
%                 than gamma in (0,1]. Default: 1. Note: changing this parameter 
%                 is not recommended.
% 'reduct'      - Boolean indicator for activating the rank-reduction
%                 technique. Default: false.
% 'threshold'   - threshold parameter for the rank-reduction technique. Two
%                 cut basis functions are flagged if the union of their active 
%                 supports relative to the union of the underlying cut elements 
%                 falls below the threshold. Default: 0.25.
%
% [Ap,P1,P2] = DEFLATION(sp_trimmed, msh_trimmed, A, int_dofs) also
% returns the left and right preconditioning matrices P1 = P and P2 = D, where
% P is the deflation preconditioner built from D*A*D and D is the diagonal
% preconditioner.
%
% [Ap,P1,P2,Pfun,Ptfun] = DEFLATION(sp_trimmed, msh_trimmed, A, int_dofs)
% also returns function handles Pfun and Ptfun for applying the projection
% and its transpose to a vector.
%
% [Ap,P1,P2,Pfun,Ptfun,r] = DEFLATION(sp_trimmed, msh_trimmed, A, int_dofs)
% also returns the deflation rank r (equal to the number of small basis
% functions among interior degrees of freedom).
%
% [Ap,P1,P2,Pfun,Ptfun,r,iss,isl] = DEFLATION(sp_trimmed, msh_trimmed, A, int_dofs)
% also returns the small and large degrees of freedom. The rows and columns
% corresponding to small degrees of freedom are zero by construction.
%
% [Ap,P1,P2,Pfun,Ptfun,r,iss,isl,Q] = DEFLATION(sp_trimmed, msh_trimmed, A, int_dofs)
% also returns the Q-operator Q(x) for recovering the approximate solution of the
% original system from the approximate solution of the deflated system and the
% right-hand side. 
%
% Reference:
% [1] J. Frank and C. Vuik. On the construction of deflation-based
% preconditioners. SISC, 2001.

%% Set algorithm parameters

Param = inputParser;
Param.KeepUnmatched = true;
Param.addRequired('sp_trimmed');
Param.addRequired('msh_trimmed');
Param.addRequired('A');
Param.addRequired('int_dofs');
Param.addParameter('gamma', 1, @(x) x>0 && x<=1);
Param.addParameter('reduct', false, @(x) isa(x, 'logical'));
Param.addParameter('threshold', 0.25, @(x) x>0 && x<=1);
Param.parse(sp_trimmed, msh_trimmed, A, int_dofs, varargin{:});

%% Retrieve parameters
gamma=Param.Results.gamma;
reduct=Param.Results.reduct;
threshold=Param.Results.threshold;

% Check if A is square
[m, n] = size(A);
if m ~= n
    error('Matrix A must be square.');
end

%% 1. Diagonal scaling
[~,D]=jacobi(A);

%% 2. Deflation based preconditioner

if reduct
    % Small basis functions that are O(1) on the trimmed elements
     IS = linear_dep(sp_trimmed, msh_trimmed, gamma, threshold);
else
    % Identify small and large basis functions (this is independent of the
    % rescaling).
    [~,~,~,IS] = msh_split(sp_trimmed, msh_trimmed, gamma);
end

% Map active to local dofs
ISr=ismember(int_dofs, IS);

I=speye(n);
Z=I(:,ISr);
Zd=D*Z;

r=sum(ISr);             % Deflation rank
Ac=Zd'*A*Zd;            % Coarse matrix
P=I-(A*Zd)*(Ac\(Zd'));  % Preconditioner (Projection)
Ap=D*P*A*D;             % Preconditioned matrix

P1=P;
P2=D;

Pfun=@(x) x-A*(Zd*(Ac\(Zd'*x)));
Ptfun=@(x) x-Zd*(Ac\(Zd'*(A*x)));

% Q-operator
Q=@(x) Zd*(Ac\(Zd'*x));

iss=find(ISr);
isl=setdiff(1:n, iss);
end