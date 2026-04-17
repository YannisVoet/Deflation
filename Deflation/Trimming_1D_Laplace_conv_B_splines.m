% Trimmed line [0, 0.75+eps]
% Conditioning of the system matrices and convergence of iterative solvers 
% for the B-spline basis for various preconditioning strategies.

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

deltarange = logspace(-4, -2.5, 20);

mrange = [7]; % refinement level
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


for j = 1:nd
    d = drange(j);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Degree %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', d);


    for i = 1:ne

        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', i);

        delta=deltarange(i);
        m=mrange;
        n = 3;
        problem_data.f = @(x) n^2 * pi^2 * sin (n*pi*x);
        problem_data.hfun = @(x, ind) sin (n*pi*x);
        problem_data.g = @(x, ind) n * pi * cos (n*pi*x);
        problem_data.uex     = @(x) sin (n*pi*x);
        problem_data.graduex = @(x) n * pi * cos (n*pi*x);

        h = 1 / (2.^m);

        method_data.nsub       = 2.^(m);
        method_data.degree     = d; % Degree of the splines
        method_data.regularity = d-1;
        % method_data.regularity = 0;
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
        rhs = op_f_v_trimming (sp_trimmed, msh_trimmed, f);

        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp_trimmed, msh_trimmed, hfun, drchlt_sides);
        int_dofs = setdiff (1:sp_trimmed.ndof, drchlt_dofs);
        rhs(int_dofs) = rhs(int_dofs) - K(int_dofs, drchlt_dofs)*u_drchlt;
        ndof=sp_trimmed.ndof;

        % Save space and mesh
        msh_cell{i}=msh_trimmed;
        sp_cell{i}=sp_trimmed;

        % Pure Neumann case
        int_dofs=1:ndof;
        n=length(int_dofs);

        Kr = K(int_dofs, int_dofs);
        Mr = M(int_dofs, int_dofs);

        [Kr,Mr]=symmetrize(Kr,Mr);

        sK=1;
        sM=0;

        %% CG parameters

        x=ones(n,1);
        b=Mr*x;
        seed=12; % Fix the random seed
        rng(seed);
        x0=rand(n,1); % Initial starting vector
        r0=Mr*x0-b;   % Initial residual

        % tol=(1e-7)*norm(b)*lambda_min; % Tolerance
        tol=1e-7; % Tolerance
        maxiter=500; % Maximum number of iterations
        param_cg.xt=x;

        %% Preconditioning strategies

        [Md,DM]=jacobi(Mr);

        for k=1:n_prec
            switch precond_tech{k}
                case 'No preconditioning'
                    [T(k).cond(i,j),lambda_min]=itcond(Mr);
                    [xt,T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Mr,b,tol*norm(b)*sqrt(lambda_min),maxiter,[],x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm(b);

                case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

                    [Md,DM]=jacobi(Mr);
                    [T(k).cond(i,j),lambda_min]=itcond(Md);
                    norm_b=normA(DM*DM,b);
                    [xt,T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Mr,b,tol*norm_b*sqrt(lambda_min),maxiter,@(x) DM*(DM*x),x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'SIPIC' % SIPIC preconditioner

                    [Ms,Ps,P2s,P1s]=sipic(Mr);
                    [T(k).cond(i,j),lambda_min]=itcond(Ms);
                    norm_b=normA(Ps*Ps',b);
                    [xt,T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Mr,b,tol*norm_b*sqrt(lambda_min),maxiter, @(x) Ps*(Ps'*x), x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'Schwarz' % Schwarz preconditioner

                    [Mz,P1z,P2z,Pfun]=schwarz(sp_trimmed, msh_trimmed, Mr, int_dofs);
                    [T(k).cond(i,j),lambda_min]=itcond(Mz);
                    norm_b=normA(@(x) DM*(Pfun(DM*x)),b);
                    [xt,T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Mr,b,tol*norm_b*sqrt(lambda_min),maxiter,@(x) DM*(Pfun(DM*x)),x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'Deflation' % Deflation based preconditioner

                    [Mp,P1p,P2p,Pfun,Ptfun,T(k).rank(i,j),~,ilM,Q]=deflation(sp_trimmed, msh_trimmed, Mr, int_dofs, T(k).param);
                    [T(k).cond(i,j),lambda_min]=itcond(Mp(ilM,ilM));
                    norm_b=normA(DM*DM,Pfun(b));
                    [xt,T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = dpcg(Mr,b,tol*norm_b*sqrt(lambda_min),maxiter,Pfun,@(x) DM*(DM*x),x0,@(x) Ptfun(x)+Q(b),param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;


                otherwise
                    error('Not implemented')
            end
            T(k).err{i,j}=T(k).err{i,j}/normA(Mr,x);
            T(k).res{i,j}=T(k).res{i,j}/norm(b);
        end


    end

end

%% Condition number

h=1/nsub;
eta_vec = deltarange/h;

% Condition number of the (preconditioned) mass
figure
for k=1:n_prec
    loglog(eta_vec, T(k).cond, 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
trendline(eta_vec, T(id_def).cond, color=T(id_def).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
title('Conditioning - $M$')
legend('Location', 'northeast')
legend show
xlabel('$\eta$')
ylabel('$\kappa$')

%% Convergence history for selected cut configuration
s=1;

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

% Error
figure
for k=1:n_prec
    semilogy(T(k).err{s,1}, 'Linestyle', '-', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Relative error')
legend('Location', 'northeast')
legend show
xlabel('Iteration number')
ylabel('Relative error')

%% Iteration numbers

figure
for k=1:n_prec
    semilogx(eta_vec, T(k).iter, 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Iteration counts - $M$')
legend('Location', 'northeast')
legend show
xlabel('$\eta$')
ylabel('Iteration count')

%% Check eigenvalues

s=1;
M = op_u_v_trimming(sp_cell{s}, sp_cell{s}, msh_cell{s});
Md=jacobi(M);
[Mz,~,~,~]=schwarz(sp_cell{s}, msh_cell{s}, M, 1:sp_cell{s}.ndof);
[Mp,~,~,~,~,rm,iss,isl]=deflation(sp_cell{s}, msh_cell{s}, M, 1:sp_cell{s}.ndof, gamma=1);

e=real(sort(eig(full(Md))));
ez=real(sort(eig(full(Mz))));
er=real(sort(eig(full((Mp(isl,isl))))));
er=[zeros(length(iss),1); er];

figure
ne=min([100 size(M,1)]);
semilogy(e(1:ne), 'Linestyle', 'none', 'LineWidth', 1, 'Marker', T(id_jac).marker, 'Color', T(id_jac).color)
grid on; hold on
semilogy(ez(1:ne), 'Linestyle', 'none', 'LineWidth', 1, 'Marker', T(id_sch).marker, 'Color', T(id_sch).color)
semilogy(er(1:ne), 'Linestyle', 'none', 'LineWidth', 1, 'Marker', T(end).marker, 'Color',  T(end).color)
title('Spectrum of $M$')
legend('Jacobi', 'Schwarz', 'Deflation', 'Location', 'Southeast')
xlabel('Eigenvalue number')
ylabel('Eigenvalue')
