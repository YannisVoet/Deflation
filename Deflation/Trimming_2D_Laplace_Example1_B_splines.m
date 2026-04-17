% Experiment 1: L^2 projection problem on a rotated lattice structure.
% Discretization with smooth B-spline bases of degree 2 or 3.
% Comparison of the deflation-based preconditioner with and without
% rank-reduction.

addpath(genpath('./Precond'))
addpath(genpath('./Auxiliary functions'))
addpath(genpath('./Geometry'))

clc
close all
clear variables

% Set the parameters
% polynomial degree
drange = [2];
% refinement level
refinement = 70;
n_ref = 0;
% case study
case_study = 'rotating_lattice_square';
eps_vec = linspace(0, pi/4, 20);

ne=length(eps_vec);
nd=length(drange);

% preconditioning techniques
precond_tech = {'No preconditioning', 'Jacobi', 'SIPIC', 'Schwarz', 'Deflation', 'Deflation'};
n_prec=length(precond_tech);
% Rank reduction
reduct=[false true];

reduct_val=num2cell(reduct);
reduct_cell=cellfun(@(x) [' (rank-reduction = ' num2str(x) ')'], reduct_val, 'UniformOutput', false);

defl_param=struct('reduct', reduct_val, 'string', reduct_cell);

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

%% Schwarz parameters

param_sch.block_sel='cut_elem'; % One index block for each cut element
param_sch.inv='exact'; % Approximate inverse for improved stability
param_sch.chol=false;

T(id_sch).param=param_sch;

%% Deflation parameters
% Rank reduction
param_defl=struct('reduct', {true, false}, 'threshold', {0.25 1});

T(end-1).param=param_defl(1);
T(end).param=param_defl(2);

T(end-1).addendum=' (reduced)';
T(end-1).color=[0.9 0.6 0];

for j = 1:length(drange)

    d = drange(j);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Degree %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', d);

    [output, square_deg] = feval(case_study,eps_vec,d,refinement);

    problem_data.nmnn_sides = [];
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

        sK=1;
        sM=0;

        %% CG parameters
        [~, ~, ~, ISl] = msh_split(sp_trimmed, msh_trimmed);

        m=4;
        uex=@(x,y) sin(m*pi*x).*sin(m*pi*y);
        b=op_f_v_trimming(sp_trimmed,msh_trimmed,uex);
        b=b(int_dofs);
        x=Mr\b;

        seed=12; % Fix the random seed
        rng(seed);
        x0=rand(n,1); % Initial starting vector
        r0=Mr*x0-b;   % Initial residual
        tol=(1e-9);   % Tolerance
        maxiter=1000; % Maximum number of iterations

        param_cg.store=false; % store intermediate solutions
        param_cg.xt=x;

        %% Preconditioning strategies

        [Md,DM]=jacobi(Mr);

        for k=1:n_prec
            switch precond_tech{k}
                case 'No preconditioning'
                    [T(k).cond(i,j),lambda_min]=itcond(Mr);
                    [T(k).sol{i,j},T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Mr,b,tol*norm(b)*sqrt(lambda_min),maxiter,[],x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm(b);

                case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

                    [Md,DM]=jacobi(Mr);
                    [T(k).cond(i,j),lambda_min]=itcond(Md);
                    norm_b=normA(DM*DM,b);
                    [T(k).sol{i,j},T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Mr,b,tol*norm_b*sqrt(lambda_min),maxiter,@(x) DM*(DM*x),x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'SIPIC' % SIPIC preconditioner

                    [Ms,Ps,P2s,P1s]=sipic(Mr);
                    [T(k).cond(i,j),lambda_min]=itcond(Ms);
                    norm_b=normA(Ps*Ps',b);
                    [T(k).sol{i,j},T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Mr,b,tol*norm_b*sqrt(lambda_min),maxiter, @(x) Ps*(Ps'*x), x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'Schwarz' % Schwarz preconditioner

                    [Mz,P1z,P2z,Pfun]=schwarz(sp_trimmed, msh_trimmed, Mr, int_dofs, T(k).param);
                    [T(k).cond(i,j),lambda_min]=itcond(Mz);
                    norm_b=normA(@(x) DM*(Pfun(DM*x)),b);
                    [T(k).sol{i,j},T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Mr,b,tol*norm_b*sqrt(lambda_min),maxiter,@(x) DM*(Pfun(DM*x)),x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'Deflation' % Deflation based preconditioner

                    [Mp,P1p,P2p,Pfun,Ptfun,T(k).rank(i,j),~,ilM,Q]=deflation(sp_trimmed, msh_trimmed, Mr, int_dofs, T(k).param);
                    [T(k).cond(i,j),lambda_min]=itcond(Mp(ilM,ilM));
                    norm_b=normA(DM*DM,Pfun(b));
                    [T(k).sol{i,j},T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = dpcg(Mr,b,tol*norm_b*sqrt(lambda_min),maxiter,Pfun,@(x) DM*(DM*x),x0,@(x) Ptfun(x)+Q(b),param_cg);
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
s=11;

% Plot trimmed configuration
msh_plot_trim_geo(output{s})

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

%% All residuals

% Non-preconditioned residual
figure
for i=1:ne
    for k=1:n_prec
        h(k)=semilogy(T(k).res{i,1}, 'Linestyle', '-', 'Marker', 'none', 'Color',  [T(k).color 0.5], 'LineWidth', 1);
        hold on; grid on;
    end
end
title('Relative residual')
legend(h, strcat({T.name}, {T.addendum}), 'Location', 'southwest')
legend show
xlim([1 400])
ylim([1e-15 1])
xlabel('Iteration number')
ylabel('Relative residual')

% Preconditioned residual
figure
for i=1:ne
    for k=1:n_prec
        h(k)=semilogy(T(k).res_prec{i,1}, 'Linestyle', '-', 'Marker', 'none', 'Color',  [T(k).color 0.5], 'LineWidth', 1);
        hold on; grid on;
    end
end
title('Relative preconditioned residual')
legend(h, strcat({T.name}, {T.addendum}), 'Location', 'southwest')
legend show
xlim([1 400])
ylim([1e-15 1])
xlabel('Iteration number')
ylabel('Relative preconditioned residual')

% Error
figure
for i=1:ne
    for k=1:n_prec
        h(k)=semilogy(T(k).err{i,1}, 'Linestyle', '-', 'Marker', 'none', 'Color',  [T(k).color 0.5], 'LineWidth', 1);
        hold on; grid on;
    end
end
title('Relative error')
legend(h, strcat({T.name}, {T.addendum}), 'Location', 'southwest')
legend show
xlim([1 400])
ylim([1e-15 1])
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

s=9;
Mr = op_u_v_trimming(sp_cell{s}, sp_cell{s}, msh_cell{s});

for k=2:n_prec

    switch precond_tech{k}
        case 'No preconditioning'

        case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

            Md=jacobi(Mr);
            T(k).e=real(sort(eig(full(Md))));

        case 'SIPIC' % SIPIC preconditioner

            Ms=sipic(Mr);
            T(k).e=real(sort(eig(full(Ms))));

        case 'Schwarz' % Schwarz preconditioner

            Mz=schwarz(sp_cell{s}, msh_cell{s}, Mr, 1:sp_cell{s}.ndof);
            T(k).e=real(sort(eig(full(Mz))));

        case 'Deflation' % Deflation based preconditioner

            [Mp,~,~,~,~,~,iss,isl]=deflation(sp_cell{s}, msh_cell{s}, Mr, 1:sp_cell{s}.ndof, T(k).param);

            er=real(sort(eig(full((Mp(isl,isl))))));
            T(k).e=[zeros(length(iss),1); er];

        otherwise
            error('Not implemented')

    end
end

nm=min([300 size(M,1)]);

figure
for k=2:n_prec
    semilogy(T(k).e(1:nm), 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Spectrum of $M$')
legend('Location', 'Southeast')
xlabel('Eigenvalue number')
ylabel('Eigenvalue')

%% Spectrum

fig=figure;
fig.Units='normalized';
fig.Position=[0 0 0.2188 0.2917];
tiledlayout(n_prec-1,1)
xmin=arrayfun(@(x) min(x.e(x.e>0)), T, 'UniformOutput', false);
xmax=arrayfun(@(x) max(x.e(x.e>0)), T, 'UniformOutput', false);

xmin=[xmin{:}]; xmin=0.5*min(xmin);
xmax=[xmax{:}]; xmax=1.5*max(xmax);

for k=2:n_prec
    nexttile
    h(k-1)=semilogx(T(k).e, zeros(1,size(Mr,1)), 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum]);
    xlim([xmin xmax])
    grid on;
end

leg=legend(h, strcat({T(2:end).name}, {T(2:end).addendum}), 'Location', 'southoutside', 'NumColumns', n_prec-1);
leg.Layout.Tile = 'South';
sgtitle('Spectrum of $M$', 'Interpreter', 'Latex', 'FontSize', 14)
fontsize(fig, 14, 'points')


%% Rank values

ranks=cat(2,T(id_defl).rank);

s=9;
msh_plot_trim_geo(output{s})

[~,~,~,IS] = msh_split(sp_cell{s}, msh_cell{s}, 1);
ISr = linear_dep(sp_cell{s}, msh_cell{s}, 1, 0.25);
msh_plot_supp(output{s}, sp_cell{s}, IS, 'red', true)
msh_plot_supp(output{s}, sp_cell{s}, ISr, 'blue', true)

%% Plot solution

npts_per_el=[5,5];
[c, data_points, connectivity]=sp_eval(T(end-1).sol{s,1}(:,end), sp_cell{s}, msh_cell{s}, geometry, npts_per_el);

xp=data_points(:,1);
yp=data_points(:,2);

cex=uex(xp,yp);

Max=max([cex c], [], 'all');
Min=min([cex c], [], 'all');

%% Visualization of solution
clear F Q Sj
% Figure appearance
F.FigParam.Units='normalized';
F.FigParam.Position=[0 0 0.4 0.3];

Q.type='patch';
Q.PlotParam=struct('Vertices', data_points, 'Faces', connectivity);
Q.AxisParam.CLim=[Min Max];
Q.AxisParam.Title.String='Exact solution';
Q.DynamicParam=struct('CData', num2cell(cex, 1));

Sj=Q;
Sj.DynamicParam=struct('CData', num2cell(c, 1));
Sj.AxisParam.Title.String='Numerical solution (Deflation)';
Sj.ColorbarParam.Visible='on';

plotsol(F,Q,Sj)

%% Exact solution

Maxex=max(cex, [], 'all');
Minex=min(cex, [], 'all');

clear F Q
% Figure appearance
F.FigParam.Units='normalized';
F.FigParam.Position=[0 0 0.2187 0.2917];
F.TiledLayout.TileSpacing='loose';
F.TiledLayout.Padding='loose';

Q.type='patch';
Q.PlotParam=struct('Vertices', data_points, 'Faces', connectivity);
Q.AxisParam.CLim=[Minex Maxex];
Q.AxisParam.XLim=[0 1];
Q.AxisParam.YLim=[0 1];
Q.AxisParam.DataAspectRatio=[1 1 1];
Q.AxisParam.Title.String='Exact solution';
Q.DynamicParam=struct('CData', num2cell(cex, 1));

plotsol(F,Q)

%% Fictitious and physical domains

s=9;
param.surf.FaceColor=[0.9 0.9 0.9];
param.surf.EdgeColor=[0.9 0.9 0.9];

figure
srf = output{s}.trim_srfs.srf;
plotnurbs(srf, [], param); % Background grid
hold on; view(2)
trimmed_srfs_plot(output{s}, 'color', [0.6 0.6 0.6], 'knots', false) % Physical domain
title('Fictitious and physical domain')
xlim([0 1])
ylim([0 1])