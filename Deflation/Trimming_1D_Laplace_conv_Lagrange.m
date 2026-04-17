% Trimmed line [0, 0.75+eps]
% Conditioning of the system matrices and convergence of iterative solvers 
% for the Lagrange basis for various preconditioning strategies.

addpath(genpath('./Precond'))
addpath(genpath('./Auxiliary functions'))
addpath(genpath('./Geometry'))

clc
close all
clear variables

% line geometry
lg = 1;
line = nrbline([0,0], [lg, 0]);
problem_data.geo_name = line;

% 1) PROBLEM DATA
problem_data.drchlt_sides = [1];
problem_data.nmnn_sides = [2];
problem_data.weak_drchlt_sides = [];

deltarange = [1e-6];

mrange = [8]; % refinement level
drange = [2]; % degree

ne=length(deltarange);
nd=length(drange);

% preconditioning techniques
precond_tech = {'No preconditioning', 'Jacobi', 'SIPIC', 'Schwarz', 'Deflation'};
n_prec=length(precond_tech);

% Load default parameters
T=def_param(precond_tech);

for k=1:n_prec
    T(k).cond=zeros(ne, nd);
    T(k).iter=zeros(ne, nd);
    T(k).res=cell(ne, nd);
    T(k).res_prec=cell(ne, nd);
    T(k).err=cell(ne, nd);
    T(k).rank=zeros(ne, nd);
end

id_def=cellfun(@(x) isequal(x,'No preconditioning'), {T.name}, 'UniformOutput', true);
id_jac=cellfun(@(x) isequal(x,'Jacobi'), {T.name}, 'UniformOutput', true);
id_sip=cellfun(@(x) isequal(x,'SIPIC'), {T.name}, 'UniformOutput', true);
id_sch=cellfun(@(x) isequal(x,'Schwarz'), {T.name}, 'UniformOutput', true);
id_defl=cellfun(@(x) isequal(x,'Deflation'), {T.name}, 'UniformOutput', true);

% Meshes and spaces
msh_cell=cell(ne,1);
sp_cell=cell(ne,1);

%% Manufactured solution
syms x

delta=deltarange;
xr=0.75+delta;
n=3;
C=8;
a=8;
w=15;
xl=1/w+xr;
r=@(x) C.^(x.^a/(xr^a)).*x.*sin(pi*1./(xl-x));

u=C.^(x.^a/(xr^a)).*x.*sin(pi*1./(xl-x));

laplacian = matlabFunction(diff(u,x,2), 'vars', [x]); % Laplacian
du_dx = matlabFunction(diff(u,x,1), 'vars', [x]);
u = matlabFunction(u, 'vars', [x]);

% Source and boundary terms
% Dimension (for vector unknown)
dim=1;

% Right-hand side
Data.f=@(x) -laplacian(x);

% Dirichlet boundary conditions
Data.g = @(x,ind) zeros ([dim, size(x)]);

% Neumann boundary conditions
Data.h=@(x,iside) du_dx(x);

% Exact solution (optional)
problem_data.uex     = u;
problem_data.graduex = du_dx;


for j = 1:nd
    d = drange(j);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Degree %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', d);


    for i = 1:ne

        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', i);

        delta=deltarange(i);
        m=mrange;
        n = 3;

        h = 1 / (2.^m);

        method_data.nsub       = 2.^(m);
        method_data.degree     = d; % Degree of the splines
        method_data.regularity = 0;
        method_data.nquad      = d+1;

        method_data.solver      = 'jacobi';

        % Manually creating the reparam structure
        if( m == 1 ) % for now handling special case separately
            warning('You should start the simulation with at least 4 elements!')
            non_trimmed_elem_ids_cell{j} = [1];
            trimmed_ids_cell{j} = [2];
        else
            non_trimmed_elem_ids_cell{j} = [1:(2^m*3/4)];
            trimmed_ids_cell{j} = [2^m*3/4 + 1];
        end

        non_trim_elem_ids = non_trimmed_elem_ids_cell{j};
        trim_elem_ids = trimmed_ids_cell{j};

        nb_non_trim_elems=length(non_trim_elem_ids);
        nb_trim_elems=length(trim_elem_ids);

        % Creating the trim_elems structure
        for k = 1:numel(trim_elem_ids)
            trim_elems(k).elem_id = non_trim_elem_ids(numel(non_trim_elem_ids))+k;
            if( m == 1 )
                % trim_elems(j).tiles = nrbsquare([0+(j-1)*element_size,0.5*length],element_size,0.25*length+delta,d,0);
            else
                reparametrization = nrbline([0.75*lg,0],[0.75*lg + delta,0]);
                trim_elems(k).tiles = nrbdegelev(reparametrization, d-1);
                trim_elems(k).nb_tiles=1;
                trim_elems(k).nb_pts=0;
                trim_elems(k).quad_pts=[];
                trim_elems(k).quad_weights=[];
            end
        end

        boundaries(1).param_side=1;
        boundaries(2).param_side=0;

        boundaries(1).label=1;
        boundaries(2).label=2;

        boundaries(1).interface=0;
        boundaries(2).interface=0;

        boundaries(1).nb_non_trim_elems=1;
        boundaries(2).nb_non_trim_elems=0;

        boundaries(1).nb_reparam_elems=0;
        boundaries(2).nb_reparam_elems=1;

        boundaries(1).non_trim_elem_bd_ids=1;
        boundaries(2).non_trim_elem_bd_ids=[];

        boundaries(1).reparam_elem_bd_ids=[];
        boundaries(2).reparam_elem_bd_ids=2.^m*3/4 + 1;

        reparam_bd=trim_elems(1);
        reparam_bd.normals=[];

        boundaries(1).reparam_elems=[];
        boundaries(2).reparam_elems=reparam_bd;

        reparam.label=-1;
        reparam.nb_non_trim_elems=nb_non_trim_elems;
        reparam.nb_trim_elems=nb_trim_elems;
        reparam.non_trim_elem_ids = non_trim_elem_ids;
        reparam.trim_elem_ids = trim_elem_ids;
        reparam.trim_elems = trim_elems;
        reparam.boundaries = boundaries;

        % Boundaries in the parametric domain
        reparam.domain=[0, 0.75 + delta];

        %% Solver

        % Extract the fields from the data structures into local variables
        data_names = fieldnames (problem_data);
        for iopt  = 1:numel (data_names)
            eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
        end
        data_names = fieldnames (method_data);
        for iopt  = 1:numel (data_names)
            eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
        end
        if ~isfield(problem_data, 'weak_drchlt_sides')
            weak_drchlt_sides = [];
        end

        % Construct geometry structure
        geometry  = geo_load(geo_name);
        [knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);

        % Construct msh structure
        rule     = msh_gauss_nodes (nquad);
        [qn, qw] = msh_set_quad_nodes (zeta, rule);
        msh_trimmed = msh_trimming (zeta, qn, qw, geometry, reparam);

        % Construct space structure
        sp_trimmed = sp_trimming (knots, degree, msh_trimmed);

        % Assemble the matrices
        K = op_gradu_gradv_trimming (sp_trimmed, sp_trimmed, msh_trimmed);
        M = op_u_v_trimming (sp_trimmed, sp_trimmed, msh_trimmed);

        % Save space and mesh
        msh_cell{i}=msh_trimmed;
        sp_cell{i}=sp_trimmed;

        % Pure Neumann boundary conditions
        [Kb,Mb]=symmetrize(K,M);

        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp_trimmed, msh_trimmed, Data.h, drchlt_sides);
        int_dofs = setdiff (1:sp_trimmed.ndof, drchlt_dofs);
        ndof=sp_trimmed.ndof;
        n=length(int_dofs);

        % Expand matrices
        ndof_full=sp_trimmed.space_untrimmed.ndof;
        active_dofs=sp_trimmed.active_dofs;

        % Change basis from Bernstein to Lagrange polynomials on
        % uniform interpolation points.
        [B,C]=interpbasis(sp_trimmed.space_untrimmed, interp_pts='uniform', trunc=false);

        C=C(active_dofs, active_dofs);
        B=B(active_dofs, active_dofs);

        K=C'*Kb*C;
        M=C'*Mb*C;

        Kr = K(int_dofs, int_dofs);
        Mr = M(int_dofs, int_dofs);

        [Kr,Mr]=symmetrize(Kr,Mr);

        %% CG parameters

        % Right-hand side vector
        f=op_f_v_gen(sp_trimmed, msh_trimmed, Data.f);
        h=op_f_v_bd_gen(sp_trimmed, msh_trimmed, nmnn_sides, Data.h);
        b=f+h;
        b=C'*b; % Change to Lagrange basis
        b=b(int_dofs);

        seed=12; % Fix the random seed
        rng(seed);
        x0=rand(n,1); % Initial starting vector
        % r0=Mr*x0-b;          % Initial residual
        % tol=(1e-7)*norm(r0); % Tolerance
        tol=(1e-7);  % Tolerance
        maxiter=1000; % Maximum number of iterations


        Ar=Kr;
        param_cg.store=true; % store intermediate solutions

        %% Preconditioning strategies

        [Ad,DA]=jacobi(Ar);

        for k=1:n_prec
            switch precond_tech{k}
                case 'No preconditioning'
                    [T(k).cond(i,j),lambda_min]=itcond(Ar);
                    [T(k).sol,T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Ar,b,tol*norm(b)*sqrt(lambda_min),maxiter,[],x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm(b);

                case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

                    [Ad,DA]=jacobi(Ar);
                    [T(k).cond(i,j),lambda_min]=itcond(Ad);
                    norm_b=normA(DA*DA,b);
                    [T(k).sol,T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Ar,b,tol*norm_b*sqrt(lambda_min),maxiter,@(x) DA*(DA*x),x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'SIPIC' % SIPIC preconditioner

                    [As,Ps,P2s,P1s]=sipic(Ar);
                    [T(k).cond(i,j),lambda_min]=itcond(As);
                    norm_b=normA(Ps*Ps',b);
                    [T(k).sol,T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Ar,b,tol*norm_b*sqrt(lambda_min),maxiter, @(x) Ps*(Ps'*x), x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'Schwarz' % Schwarz preconditioner

                    [Az,P1z,P2z,Pfun]=schwarz(sp_trimmed, msh_trimmed, Ar, int_dofs, T(k).param);
                    [T(k).cond(i,j),lambda_min]=itcond(Az);
                    norm_b=normA(@(x) DA*(Pfun(DA*x)),b);
                    [T(k).sol,T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Ar,b,tol*norm_b*sqrt(lambda_min),maxiter,@(x) DA*(Pfun(DA*x)),x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'Deflation' % Deflation based preconditioner

                    [Ap,P1p,P2p,Pfun,Ptfun,T(k).rank(i,j),~,ilM,Q]=deflation(sp_trimmed, msh_trimmed, Ar, int_dofs, T(k).param);
                    [T(k).cond(i,j),lambda_min]=itcond(Ap(ilM,ilM));
                    norm_b=normA(DA*DA,Pfun(b));
                    [T(k).sol,T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = dpcg(Ar,b,tol*norm_b*sqrt(lambda_min),maxiter,Pfun,@(x) DA*(DA*x),x0,@(x) Ptfun(x)+Q(b),param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;


                otherwise
                    error('Not implemented')
            end
            T(k).res{i,j}=T(k).res{i,j}/norm(b);
        end


    end

end

%% Post-processing
s=1;

rhs=op_f_v_gen(sp_trimmed, msh_trimmed, problem_data.uex);
M = op_u_v_trimming(sp_cell{s}, sp_cell{s}, msh_cell{s});
v_ref=M\rhs;

for k=1:n_prec
    v=zeros(ndof,T(k).iter+1);
    v(int_dofs,:)=T(k).sol;
    % Back-transform to the B-spline basis
    v=C*v;
    T(k).sol=v;
    T(k).errl2{s,1}=sp_l2_error(sp_trimmed, msh_trimmed, T(k).sol, uex);
    T(k).errl2_proj{s,1}=sqrt(sum((T(k).sol-v_ref).*(M*(T(k).sol-v_ref)),1));
end

best_errl2=sp_l2_error(sp_trimmed, msh_trimmed, v_ref, uex);

%% Condition number

% h=1/nsub;
% eta_vec = deltarange/h;
% 
% % Condition number of the (preconditioned) mass
% figure
% for k=1:n_prec
%     loglog(eta_vec, T(k).cond, 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
%     hold on; grid on;
% end
% trendline(eta_vec, T(id_def).cond, color=T(id_def).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
% title('Conditioning - $M$')
% legend('Location', 'northeast')
% legend show
% xlabel('$\eta$')
% ylabel('$\kappa$')

%% Convergence history for selected cut configuration

% Non-preconditioned residual
figure
for k=1:n_prec
    semilogy(T(k).res{s,1}, 'Linestyle', '-', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Relative residual')
legend('Location', 'northeast')
legend show
xlabel('Iteration number')
ylabel('Relative residual')

% Preconditioned residual
figure
for k=1:n_prec
    semilogy(T(k).res_prec{s,1}, 'Linestyle', '-', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Relative preconditioned residual')
legend('Location', 'northeast')
legend show
xlabel('Iteration number')
ylabel('Relative preconditioned residual')

% L^2 error
figure
for k=1:n_prec
    semilogy(T(k).errl2{s,1}, 'Linestyle', '-', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('$L^2$ Error')
legend('Location', 'northeast')
legend show
xlabel('Iteration number')
ylabel('$L^2$ Error')

% L^2 projection error
figure
for k=1:n_prec
    semilogy(T(k).errl2_proj{s,1}, 'Linestyle', '-', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('$L^2$ projection error')
legend('Location', 'northeast')
legend show
xlabel('Iteration number')
ylabel('$L^2$ projection error')

%% Iteration numbers

% figure
% for k=1:n_prec
%     semilogx(eta_vec, T(k).iter, 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
%     hold on; grid on;
% end
% title('Iteration counts - $M$')
% legend('Location', 'northeast')
% legend show
% xlabel('$\eta$')
% ylabel('Iteration count')

%% Check eigenvalues

s=1;
K = op_gradu_gradv_trimming(sp_cell{s}, sp_cell{s}, msh_cell{s});
K=C'*K*C;

Ar=K(int_dofs,int_dofs);

for k=2:n_prec

    switch precond_tech{k}
        case 'No preconditioning'

        case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

            Ad=jacobi(Ar);
            T(k).e=real(sort(eig(full(Ad))));

        case 'SIPIC' % SIPIC preconditioner

            As=sipic(Ar);
            T(k).e=real(sort(eig(full(As))));

        case 'Schwarz' % Schwarz preconditioner

            Az=schwarz(sp_cell{s}, msh_cell{s}, Ar, int_dofs);
            T(k).e=real(sort(eig(full(Az))));

        case 'Deflation' % Deflation based preconditioner

            [Ap,~,~,~,~,~,iss,isl]=deflation(sp_cell{s}, msh_cell{s}, Ar, int_dofs, T(k).param);

            er=real(sort(eig(full((Ap(isl,isl))))));
            T(k).e=[zeros(length(iss),1); er];

        otherwise
            error('Not implemented')

    end
end

ne=min([100 size(Ar,1)]);

figure
for k=2:n_prec
    semilogy(T(k).e(1:ne), 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Spectrum of $M$')
legend('Location', 'Southeast')
xlabel('Eigenvalue number')
ylabel('Eigenvalue')

%% Rank values
ranks=cat(2,T(id_defl).rank);

%% Components

VK=zeros(ndof,n);

[VK(int_dofs,:), eigK]=eig(full(symmetrize(Kr)), 'vector');

% Back-transforming the eigenvectors
VK(int_dofs,:)=d.*VK(int_dofs,:);

[eigK, IKM]=sort(eigK);

VK=VK(:,IKM);
alpha_v_VK=(VK'*M*VK)\VK'*rhs;  % Coefficients in the eigenbasis

%% Direct solve

vc=zeros(ndof,1);
vc(int_dofs)=Kr\b;
vc=C*vc;

%% Visualize final solution
npts_per_el=[20,1];
X=[T(id_jac).sol(:,end) T(id_defl).sol(:,end)];
[c, data_points, connectivity]=sp_eval(X, sp_cell{s}, msh_cell{s}, geometry, npts_per_el);

c_jac=c(:,1);
c_defl=c(:,2);

xp=data_points(:,1);
cex=uex(xp);

Max=max([cex c_jac c_defl], [], 'all');
Min=min([cex c_jac c_defl], [], 'all');

dp=cat(3, cex, c_jac, c_defl);
dp=permute(dp, [1 3 2]);
dp=num2cell(dp, [1 2]);
dp=squeeze(dp);

%% Visualization of solution
clear F Sd
% Figure appearance
F.FigParam.Units='normalized';
F.FigParam.Position=[0 0 0.6 0.3];

% Plot all curves
Sd.type='plot';
Sd.PlotParam.XData=xp;
Sd.PlotParam.Color={'k', T(id_jac).color, T(id_defl).color};
Sd.AxisParam.XLim=[0 0.8];
Sd.AxisParam.YLim=[Min Max];
Sd.AxisParam.Title.String='Displacement';
Sd.LegendParam.Visible='on';
Sd.LegendParam.Location='southwest';
Sd.LegendParam.Title.String={'Exact solution', T(id_jac).name, T(id_defl).name};
Sd.DynamicParam=struct('YData', dp);

plotsol(F,Sd)

%% Exact solution

Maxex=max(cex, [], 'all');
Minex=min(cex, [], 'all');

clear F Sd
% Figure appearance
F.FigParam.Units='normalized';
F.FigParam.Position=[0 0 0.6 0.3];

% Plot all curves
Sd.type='plot';
Sd.PlotParam.XData=xp;
Sd.PlotParam.Color='k';
Sd.AxisParam.XLim=[0 0.8];
Sd.AxisParam.YLim=[Minex Maxex];
Sd.AxisParam.Title.String='Exact solution';
Sd.LegendParam.Visible='on';
Sd.LegendParam.Location='southwest';
Sd.LegendParam.Title.String={'Exact solution'};
Sd.DynamicParam=struct('YData', cex);

plotsol(F,Sd)
