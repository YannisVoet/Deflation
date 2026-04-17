function[Ks,varargout]=sipic(K, varargin)

% SIPIC: implements the SIPIC preconditioner from [1].
%
% Ks = SIPIC(K) computes the SIPIC preconditioned matrix Ks = G'*K*G where
% G is generally a composition of a diagonal scaling D (Jacobi
% preconditioner) and a sparse matrix S resulting from local
% orthonormalizations.
%
% Ks = SIPIC(K,gamma) specifies the orthonormalization threshold, which
% takes real positive values in [0,1]. Orthgonalization is applied if an
% off-diagonal entry of |D'*K*D| is larger than gamma. Small values of
% gamma improve the quality of the preconditioner but also increase the
% fill-in. Default: 0.9.
%
% Ks = SIPIC(K,gamma,tol) also specifies the tolerance for detecting linear
% dependencies and takes real positive values in (eps,1). An off-diagonal
% entry of |D'*K*D| between 1-tol and 1-eps indicates linear dependency
% beyond repair. This parameter allows to preemptively eliminate rows and
% columns of K that will cause numerical issues during orthonormalization.
% Default: 1e2*eps ≈ 2e-14.
%
% [Ks,G] = SIPIC(K,...) returns the preconditioning matrix G = D*S.
%
% [Ks,G,D,S] = SIPIC(K,...) also returns the diagonal scaling matrix D and 
% orthonormalization matrix S.
%
% [Ks,G,D,S,s1,s2] = SIPIC(K,...) also returns the indices s1 and s2 for 
% the nonsingular and singular submatrices of K.
%
% Reference:
% [1] F. de Prenter, C. V. Verhoosel, G. J. van Zwieten, and E. H. 
% van Brummelen. Condition number analysis and preconditioning of the 
% finite cell method. CMAME, 2017.

%% Set algorithm parameters

% Set default parameters (see [1])
Default{1}=0.9;
Default{2}=1e2*eps;

% Replace empty inputs with default parameters
def=cell2mat(cellfun(@isempty, varargin, 'UniformOutput', false));
[varargin{def}]=Default{def};

n=size(K,1);
Param = inputParser;
Param.addRequired('K', @(x) all(size(x)==[n n]));
Param.addOptional('gamma', Default{1}, @(x) isscalar(x) && x > 0);
Param.addOptional('tol', Default{2}, @(x) isscalar(x) && x > 0);
Param.parse(K, varargin{:});

%% Retrieve parameters
gamma=Param.Results.gamma;
tol=Param.Results.tol;

% First diagonal scaling
d = 1./sqrt(diag(K));

% Diagonally scaled matrix (correlation matrix)
Ks=d.*K.*d';

% Pre-processing step: remove rows and columns that are linearly dependent
% up to machine precision
[Is,Js]=find((1-tol)<abs(Ks) & abs(Ks)<1-eps);

% Form graph and identify connected components
Gs=graph(Is,Js);
s2=conncomp(Gs,'OutputForm','cell');

idx=cellfun(@length, s2);
s2=s2(idx>1);

% Eliminate rows with the smallest norm
if ~isempty(s2)
    for k=1:length(s2)
        [~,idx]=sort(vecnorm(K(s2{k},:),2,2), 'descend');
        s2{k}=s2{k}(idx(2:end));
    end
end

s2=[s2{:}];
s1=setdiff(1:n,s2);

K=K(s1,s1);
Ks=Ks(s1,s1);
d=d(s1);
n=length(s1);

% Find large off-diagonal entries
[I,J]=find(abs(Ks-speye(n,n))>gamma);

% Form graph and identify connected components
G=graph(I,J);
c=conncomp(G,'OutputForm','cell');

idx=cellfun(@length, c);
c=c(idx>1);

% Initialization
D=spdiags(d,0,n,n);
Si=speye(n,n);
Sf=speye(n,n);
niter=1;
maxiter=10;

while ~isempty(c) && niter<maxiter
    nc=length(c);

    % Sort indices according to the least nonzero row entries in K
    for k=1:nc
        [~,idx]=sort(sum(abs(K(c{k},:))>0,2));
        c{k}=c{k}(idx);

        Si(c{k},c{k})=mgs(diag(diag(K(c{k},c{k}))),Ks(c{k},c{k}));
    end

    Sf=Sf*Si;
    Ks=Sf'*D*K*D*Sf;

    % Find large off-diagonal entries
    [I,J]=find(abs(Ks-speye(n,n))>gamma);

    % Form graph and identify connected components
    G=graph(I,J);
    c=conncomp(G,'OutputForm','cell');

    idx=cellfun(@length, c);
    c=c(idx>1);
    niter=niter+1;
    Si=speye(n,n);
end

varargout{1}=D*Sf;
varargout{2}=D;
varargout{3}=Sf;
varargout{4}=s1;
varargout{5}=s2;

end