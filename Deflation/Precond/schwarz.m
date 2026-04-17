function[Ap, P1, P2, fun]=schwarz(sp_trimmed, msh_trimmed, A, int_dofs, varargin)

% SCHWARZ: Computes the additive Schwarz preconditioner defined in [1].
%
% Ap = SCHWARZ(msh_trimmed, sp_trimmed, A, int_dofs) returns the
% preconditioned matrix for the composition of diagonal scaling and Schwarz
% preconditioning. In the non-symmetric case, Ap = P*D*A*D, where 
% D = 1./sqrt(diag(A)) and P is the Schwarz preconditioner. In the symmetric 
% case, Ap = L'*D*A*D*L, where L is the lower triangular Cholesky factor of P.
%
% [Ap, P1, P2] = SCHWARZ(msh_trimmed, sp_trimmed, A, int_dofs) also returns
% the left and right preconditioners P1 = P and P2 = D, where D is the
% diagonal preconditioner for A and P is the Schwarz preconditioner for D*A*D.
%
% [Ap, P1, P2, Pfun] = SCHWARZ(msh_trimmed, sp_trimmed, A, int_dofs) also
% returns the function handle for computing P*x.
%
% [Ap, P1, P2, Pfun] = SCHWARZ(msh_trimmed, sp_trimmed, A, int_dofs, options)
% specifies optional name/value pair arguments, including:
% 'block_sel':  - specifies the block selection strategy. Available strategies are:
%       'cut_elem'    - forms an index block for each cut element by including
%                       the indices of all basis functions that are supported on
%                       that element (default).
%       'bad_func'    - forms a single index block containing the indices of 
%                       all basis functions that are only supported on cut 
%                       elements.
%       'overlap'     - forms an index block for each cut basis function
%                       containing the indices of all other functions that 
%                       intersect it.
% 'inv':        - specifies the strategy for computing the inverses of the
%                 submatrices. Available strategies are:
%       'exact'       - exact computation of inverses with inv (less stable,
%                       but more general, default).
%       'approx'      - approximate computation of inverses via their
%                       truncated spectral decomposition (more stable,
%                       but only for symmetric matrices, see [2]).
% 'tol':        - truncation tolerance for computing the approximate inverses. 
%                 All eigenvalues smaller than tol are discarded. Default: 1e-14.
% 'chol':       - if true, computes a symmetric preconditioned matrix using the 
%                 Cholesky factorization of the preconditioning matrix P. 
%                 Default: true. Note: if this parameter is false the 
%                 preconditioned matrix is non-symmetric, even if A is.
%
% Reference:
% [1] F. de Prenter, C. Verhoosel, and E. Van Brummelen. Preconditioning
% immersed isogeometric finite element methods with application to flow
% problems. CMAME, 2019.
% [2] J. N. Jomo, F. de Prenter, M. Elhaddad, D. D2Angella, C. V. Verhoosel, 
% S. Kollmannsberger, J. S. Kirschke, V. N¨ubel, E. van Brummelen, and 
% E. Rank. Robust and parallel scalable iterative solutions for large-scale 
% finite cell analyses. Finite Elements in Analysis and Design, 2019.

%% Set algorithm parameters

Param = inputParser;
Param.KeepUnmatched = true;
Param.addRequired('sp_trimmed');
Param.addRequired('msh_trimmed');
Param.addRequired('A');
Param.addRequired('int_dofs');
Param.addParameter('block_sel', 'cut_elem', @(x) ismember(x,{'cut_elem','bad_func','overlap'}));
Param.addParameter('inv', 'exact', @(x) ismember(x,{'exact','approx'}));
Param.addParameter('tol', 1e-14, @(x) x>0);
Param.addParameter('chol', true, @(x) isa(x, 'logical'));
Param.parse(sp_trimmed, msh_trimmed, A, int_dofs, varargin{:});

%% Retrieve parameters
block_sel=Param.Results.block_sel;
inverse=Param.Results.inv;
tol=Param.Results.tol;
cchol=Param.Results.chol;

%% Block selection strategy

switch block_sel
    case 'cut_elem'
        % Indices all trimmed elements
        elem_ids = msh_trimmed.reparam.trim_elem_ids;
        % Indices blocks extraction
        [~,M_vec] = sp_get_basis_functions(sp_trimmed, msh_trimmed, elem_ids);

    case 'bad_func'
        % Consider all trimmed elements
        gamma=1;
        % Indices of small basis functions
        [~,~,~,M_vec{1}] = msh_split(sp_trimmed, msh_trimmed, gamma);

    case 'overlap'
        % Consider all trimmed elements
        gamma=1;
        % Indices of small basis functions
        [~,~,~,ISl] = msh_split(sp_trimmed, msh_trimmed, gamma);

        S=op_u_v_trimming (sp_trimmed, sp_trimmed, msh_trimmed);
        S=abs(S)>0;

        nf=length(ISl);
        M_vec=cell(nf,1);

        for k=1:nf
            M_vec{k}=find(S(:,ISl(k)));
        end
end

%% 1. Diagonal preconditioning
[As,D]=jacobi(A);

%% 2. Schwarz preconditioning

bad_dofs=unique(cat(1,M_vec{:}));
good_dofs=setdiff(int_dofs, bad_dofs)';
S_vec=num2cell(good_dofs);

K_vec=[M_vec; S_vec];

% Map global indices to local ones
K_vec=cellfun(@(x) find(ismember(int_dofs', x)), K_vec, 'UniformOutput', false);
nv=length(K_vec);

A_inv=cell(nv,1);

% Loop over index blocks
for i = 1:nv

    switch inverse
        case 'exact'
            A_inv{i}=inv(As(K_vec{i},K_vec{i}));

        case 'approx'
            [Us,Ds]=eig(full(As(K_vec{i},K_vec{i})));
            d=diag(Ds);
            s=abs(d)>tol;
            A_inv{i}=(Us(:,s)./d(s)')*Us(:,s)';
    end

end

P=assemble(A_inv, K_vec);
[L, flag]=chol(P,'lower');

% Preconditioned matrix
if flag==0 && cchol % P is numerically positive definite
    Ap=L'*As*L;
else % Form a similar non-symmetric preconditioned matrix
    Ap = P*As;
end

% Preconditioning matrices
P1=P;
P2=D;

% Function handle for P*x
fun = @(x) Pfun(A_inv, K_vec, x);

    function[rhs] = Pfun(A_inv, K_vec, x)

        xb=cellfun(@(y) x(y), K_vec, 'UniformOutput', false);
        M=cellfun(@(x,y) x*y, A_inv, xb, 'UniformOutput', false);
        V=cat(1, M{:});
        I=cat(1, K_vec{:});
        rhs=accumarray(I,V);
    end

end