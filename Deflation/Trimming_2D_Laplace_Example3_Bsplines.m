% Wave equation on a waveguide-inspired spiky structure.
% Simulation of the first time step of an implicit solver.
% Performance comparison of various preconditioning strategies.

addpath(genpath('./Precond'))
addpath(genpath('./Auxiliary functions'))
addpath(genpath('./Geometry'))

clc
close all
clear variables

% Set the parameters
% polynomial degree
drange = [3];
% refinement level
refinement = [32 48];
n_ref = 0;
% case study
case_study = 'waveguide';
eps_vec = logspace(-2.5,-3.5,20);

ne=length(eps_vec);
nd=length(drange);

% preconditioning techniques
precond_tech = {'Jacobi', 'SIPIC', 'Schwarz', 'Deflation', 'Deflation'};
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

% vector with the smallest volume ratio
eta_vec = zeros(ne,1);

% Override defaults
%% Deflation parameters
% Rank reduction
param_defl=struct('reduct', {true, false}, 'threshold', {0.25 1});

T(end-1).param=param_defl(1);
T(end).param=param_defl(2);

T(end-1).addendum=' (reduced)';
T(end-1).color=[0.9 0.6 0];

%% Schwarz parameters

% param_sch.block_sel='cut_elem'; % One index block for each cut element
param_sch.block_sel='overlap'; % block index selection based on overlapping basis functions.
param_sch.inv='approx'; % Approximate inverse for improved stability

T(id_sch).param=param_sch;

%% Set data

Data.f=@(x,y,t) zeros ([1, size(x)]);

% Dirichlet boundary conditions
Data.g = @(x,y,t,ind) zeros ([1, size(x)]);
Data.d1_g = @(x,y,t,ind) zeros ([1, size(x)]);
Data.d2_g = @(x,y,t,ind) zeros ([1, size(x)]);

% Neumann boundary conditions
Data.h = @(x,y,t,ind) zeros ([1, size(x)]);

% Initial conditions
sigma=1e-1;
Data.u0 = @(x,y) exp(-((x-0.4844)/sigma).^2);
Data.d1_u0= @(x,y) zeros ([1, size(x)]);


for j = 1:length(drange)

    d = drange(j);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Degree %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', d);

    [output, square_deg] = feval(case_study,eps_vec,d,refinement);

    problem_data.nmnn_sides = [1:14];
    problem_data.drchlt_sides = [];
    problem_data.weak_drchlt_sides = [];

    problem_data.hfun = @(x, y, ind) zeros(size(x));

    method_data.stabilization = false;
    method_data.stabilization_type = 'physical';
    method_data.theta = 1e-1;
    method_data.Nitsche_type = 'symmetric';
    method_data.Cpen = (d+1)^2*10;

    method_data.nsub       = 2.^[n_ref n_ref];
    method_data.degree     = [d d];     % Degree of the splines
    method_data.regularity = [d-1 d-1]; % Regularity of the splines
    % method_data.regularity = [0 0]; % Regularity of the splines
    method_data.nquad      = [d+1 d+1];

    % Geometry definition
    problem_data.geo_name = square_deg;

    for i = 1:ne
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', i);

        method_data.reparam = output{i}(1).trim_srfs(1);

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

        % Number of subdivisions in each parametric direction
        N=cellfun(@length, zeta, 'UniformOutput', true)-1;

        % Construct msh structure
        rule     = msh_gauss_nodes (nquad);
        [qn, qw] = msh_set_quad_nodes (zeta, rule);
        msh_trimmed = msh_trimming (zeta, qn, qw, geometry, reparam);

        % Construct space structure
        sp_trimmed = sp_trimming (knots, degree, msh_trimmed);

        M = op_u_v_trimming (sp_trimmed, sp_trimmed, msh_trimmed);
        K = op_gradu_gradv_trimming (sp_trimmed, sp_trimmed, msh_trimmed);

        % Get the volume of the smallest trimmed element
        [vol_trim, size_trim] = msh_get_element_size(msh_trimmed);
        eta_vec(i) = min(vol_trim)/size_trim(1);

        % Save space and mesh
        msh_cell{i}=msh_trimmed;
        sp_cell{i}=sp_trimmed;

        % Pure Neumann boundary conditions
        [K,M]=symmetrize(K,M);

        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp_trimmed, msh_trimmed, hfun, drchlt_sides);
        int_dofs = setdiff (1:sp_trimmed.ndof, drchlt_dofs);
        ndof=sp_trimmed.ndof;
        n=length(int_dofs);

        Kr = K(int_dofs, int_dofs);
        Mr = M(int_dofs, int_dofs);

        [Kr,Mr]=symmetrize(Kr,Mr);

        %% Right-hand side for computing u1

        beta=1/4;
        gamma=1/2;
        % Time span
        I=[0 0.4];
        N=600;
        Dt=(I(2)-I(1))/N;

        F0=op_f_v_gen(sp_trimmed, msh_trimmed, Data.u0);
        U0=M\F0;

        A0=M\(-K*U0);
        Ut1=U0+(1/2-beta)*Dt^2*A0;
        Vt1=(1-gamma)*Dt*A0;
        Ar=Mr+beta*Dt^2*Kr;
        b=-K*Ut1;

        % For improved accuracy
        [Ad,DA]=jacobi(Ar);
        y=Ad\(DA*b);
        x=DA*y;

        %% CG parameters

        seed=12; % Fix the random seed
        rng(seed);
        x0=rand(n,1); % Initial starting vector
        tol=(1e-9);   % Tolerance
        maxiter=1e3;  % Maximum number of iterations

        param_cg.store=false; % store intermediate solutions
        param_cg.xt=x;

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
            T(k).err{i,j}=T(k).err{i,j}/normA(Ar,x);
            T(k).res{i,j}=T(k).res{i,j}/norm(b);
        end


    end

end

%% Condition number

% Condition number of the (preconditioned) mass
figure
for k=1:n_prec
    loglog(eta_vec, T(k).cond, 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
trendline(eta_vec(5:end), T(id_jac).cond(5:end), color=T(id_jac).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
trendline(eta_vec, T(id_sch).cond, color=T(id_sch).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
title('Conditioning - $M+\beta\Delta t^2 K$')
legend('Location', 'northeast')
legend show
xlabel('$\eta$')
ylabel('$\kappa$')

%% Convergence history for selected cut configuration
s=1;

% Plot trimmed configuration
msh_plot_trim_geo(output{s})

% Non-preconditioned residual
figure
for k=1:n_prec
    semilogy(T(k).res{s,1}, 'Linestyle', '-', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Convergence history - $M + \beta \Delta t^2 K$')
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
title('Convergence history - $M + \beta \Delta t^2 K$')
legend('Location', 'northeast')
legend show
xlabel('Iteration number')
ylabel('Preconditioned relative residual')

% Error
figure
for k=1:n_prec
    semilogy(T(k).err{s,1}, 'Linestyle', '-', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Convergence history - $M + \beta \Delta t^2 K$')
legend('Location', 'northeast')
legend show
xlabel('Iteration number')
ylabel('$A$-norm error')

%% All residuals

% Non-preconditioned residual
figure
for i=1:ne
    for k=1:n_prec
        p(k)=semilogy(T(k).res{i,1}, 'Linestyle', '-', 'Marker', 'none', 'Color',  [T(k).color 0.5], 'LineWidth', 1);
        hold on; grid on;
    end
end

title('Relative residual')
legend(p, strcat({T.name}, {T.addendum}), 'Location', 'southeast')
legend show
xlim([1 1e3])
xlabel('Iteration number')
ylabel('Relative residual')

% Preconditioned residual
figure
for i=1:ne
    for k=1:n_prec
        p(k)=semilogy(T(k).res_prec{i,1}, 'Linestyle', '-', 'Marker', 'none', 'Color',  [T(k).color 0.5], 'LineWidth', 1);
        hold on; grid on;
    end
end
title('Relative preconditioned residual')
legend(p, strcat({T.name}, {T.addendum}), 'Location', 'southeast')
legend show
xlim([1 1e3])
xlabel('Iteration number')
ylabel('Relative preconditioned residual')

% Error
figure
for i=1:ne
    for k=1:n_prec
        p(k)=semilogy(T(k).err{i,1}, 'Linestyle', '-', 'Marker', 'none', 'Color',  [T(k).color 0.5], 'LineWidth', 1);
        hold on; grid on;
    end
end
title('Relative error')
legend(p, strcat({T.name}, {T.addendum}), 'Location', 'southeast')
legend show
xlim([1 1e3])
xlabel('Iteration number')
ylabel('Relative error')

%% Iteration numbers

% System size
ndofs=cellfun(@(x) x.ndof, sp_cell, 'UniformOutput', true);

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

s=20;
K = op_gradu_gradv_trimming(sp_cell{s}, sp_cell{s}, msh_cell{s});
M = op_u_v_trimming(sp_cell{s}, sp_cell{s}, msh_cell{s});

Kr = K(int_dofs, int_dofs);
Mr = M(int_dofs, int_dofs);
Ar=Mr+beta*Dt^2*Kr;

for k=1:n_prec

    switch precond_tech{k}
        case 'No preconditioning'

        case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

            Ad=jacobi(Ar);
            T(k).e=real(sort(eig(full(Ad))));

        case 'SIPIC' % SIPIC preconditioner

            As=sipic(Ar);
            T(k).e=real(sort(eig(full(As))));

        case 'Schwarz' % Schwarz preconditioner

            Az=schwarz(sp_cell{s}, msh_cell{s}, Ar, int_dofs, T(k).param);
            T(k).e=real(sort(eig(full(Az))));

        case 'Deflation' % Deflation based preconditioner

            [Ap,~,~,~,~,~,iss,isl]=deflation(sp_cell{s}, msh_cell{s}, Ar, int_dofs, T(k).param);

            er=real(sort(eig(full((Ap(isl,isl))))));
            T(k).e=[zeros(length(iss),1); er];

        otherwise
            error('Not implemented')

    end
end

nr=min([100 size(Ar,1)]);

figure
for k=1:n_prec
    semilogy(T(k).e(1:nr), 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Spectrum of $M + \beta \Delta t^2 K$')
legend('Location', 'Southeast')
xlabel('Eigenvalue number')
ylabel('Eigenvalue')

%% Spectrum

clear p
fig=figure;
fig.Units='normalized';
fig.Position=[0 0 0.2188 0.2917];
tiledlayout(n_prec,1)
xmin=arrayfun(@(x) min(x.e(x.e>0)), T, 'UniformOutput', false);
xmax=arrayfun(@(x) max(x.e(x.e>0)), T, 'UniformOutput', false);

xmin=[xmin{:}]; xmin=0.5*min(xmin);
xmax=[xmax{:}]; xmax=1.5*max(xmax);
nk=ndof;

for k=1:n_prec
    nexttile
    p(k)=semilogx(T(k).e, zeros(1,nk), 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum]);
    xlim([xmin xmax])
    grid on;
end

leg=legend(p, strcat({T.name}, {T.addendum}), 'Location', 'southoutside', 'NumColumns', n_prec);
leg.Layout.Tile = 'South';
sgtitle('Spectrum of $M + \beta \Delta t^2 K$', 'Interpreter', 'Latex', 'FontSize', 14)
fontsize(fig, 14, 'points')

%% Rank values
ranks=cat(2,T(id_defl).rank);

s=20;
msh_plot_trim_geo(output{s})

[~,~,~,IS] = msh_split(sp_cell{s}, msh_cell{s}, 1);
ISr = linear_dep(sp_cell{s}, msh_cell{s}, 1, 0.25);
msh_plot_supp(output{s}, sp_cell{s}, IS, 'red', true)
msh_plot_supp(output{s}, sp_cell{s}, ISr, 'blue', true)
