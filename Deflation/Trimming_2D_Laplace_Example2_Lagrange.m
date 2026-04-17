% Experiment 2: Poisson problem on a extruded domain discretized with
% quadratic Lagrange bases. 

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
refinement = 56;
n_ref = 0;
% case study
case_study = 'extrusion';
eps_vec = logspace(-2,-4,20);

ne=length(eps_vec);
nd=length(drange);

% preconditioning techniques
precond_tech = {'Jacobi', 'SIPIC', 'Schwarz', 'Deflation'};
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
C_cell=cell(ne,1);
int_dofs_cell=cell(ne,1);

% vector with the smallest volume ratio
eta_vec = zeros(ne,1);

% Override defaults
%% Schwarz parameters

param_sch.block_sel='cut_elem'; % One index block for each cut element
param_sch.inv='approx'; % Approximate inverse for improved stability

T(id_sch).param=param_sch;

%% Manufactured solution

syms x y

u=x*(1-x)*(sin(3*pi*x))^2.*sin(pi*y);

laplacian = matlabFunction(diff(u,x,2)+diff(u,y,2), 'vars', [x y]); % Laplacian
du_dx = matlabFunction(diff(u,x,1), 'vars', [x y]);
du_dy = matlabFunction(diff(u,y,1), 'vars', [x y]);
u = matlabFunction(u, 'vars', [x y]);

x_fine=linspace(0,1,1000);
y_fine=linspace(0,1,1000);
[X,Y] = meshgrid(x_fine, y_fine);

figure
mesh(X,Y,u(X,Y))
colormap('jet');
view(2)

Data.f= @(x,y) -laplacian(x,y);

% Dirichlet boundary conditions
Data.g = @(x,y,ind) zeros ([1, size(x)]);

c5=[0.5,0.25];
c7=[0.5,0.75];

% Neumann boundary conditions
v51 = @(x,y) -(x-c5(1))./sqrt((x-c5(1)).^2+(y-c5(2)).^2);
v52 = @(x,y) -(y-c5(2))./sqrt((x-c5(1)).^2+(y-c5(2)).^2);

v71 = @(x,y) -(x-c7(1))./sqrt((x-c7(1)).^2+(y-c7(2)).^2);
v72 = @(x,y) -(y-c7(2))./sqrt((x-c7(1)).^2+(y-c7(2)).^2);

% Outward pointing normals
v6=[-1; 0]; v3=[0 -1];
v8=[1; 0];  v4=[0 1];

grad_u = @(x,y) [du_dx(x,y); du_dy(x,y)];
Data.h=@(x,y,iside) (du_dx(x,y)*v3(1)+du_dy(x,y)*v3(2))*(iside==3)+(du_dx(x,y)*v4(1)+du_dy(x,y)*v4(2))*(iside==4)+...
    (du_dx(x,y)*v6(1)+du_dy(x,y)*v6(2))*(iside==6)+(du_dx(x,y)*v8(1)+du_dy(x,y)*v8(2))*(iside==8)+...
    (du_dx(x,y).*v51(x,y)+du_dy(x,y).*v52(x,y))*(iside==5)+(du_dx(x,y).*v71(x,y)+du_dy(x,y).*v72(x,y))*(iside==7);

% Exact solution (optional)
problem_data.uex     = u;
problem_data.graduex = @(x,y) cat(1, ...
    reshape (du_dx (x,y), [1, size(x)]), ...
    reshape (du_dy (x,y), [1, size(x)]));

for j = 1:length(drange)

    d = drange(j);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Degree %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', d);

    [output, square_deg] = feval(case_study,eps_vec,d,refinement);

    problem_data.nmnn_sides = [3 4 5 6 7 8];
    problem_data.drchlt_sides = [1 2];
    problem_data.weak_drchlt_sides = [];

    problem_data.hfun = @(x, y, ind) zeros(size(x));

    method_data.stabilization = false;
    method_data.stabilization_type = 'physical';
    method_data.theta = 1e-1;
    method_data.Nitsche_type = 'symmetric';
    method_data.Cpen = (d+1)^2*10;

    method_data.nsub       = 2.^[n_ref n_ref];
    method_data.degree     = [d d]; % Degree of the splines
    method_data.regularity = [0 0]; % Regularity of the splines
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
        [Kb,Mb]=symmetrize(K,M);

        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp_trimmed, msh_trimmed, hfun, drchlt_sides);
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

        C_cell{i}=C;
        int_dofs_cell{i}=int_dofs;

        %% CG parameters

        % Right-hand side vector
        f=op_f_v_gen(sp_trimmed, msh_trimmed, Data.f);
        h=op_f_v_bd_gen(sp_trimmed, msh_trimmed, nmnn_sides, Data.h);
        b=f+h;
        b=C'*b; % Change to Lagrange basis
        b=b(int_dofs);

        seed=12; % Fix the random seed
        rng(seed);
        x0=rand(n,1);  % Initial starting vector
        x=Kr\b;
        tol=(1e-9);    % Tolerance
        maxiter=10000; % Maximum number of iterations

        Ar=Kr;
        param_cg.store=true; % store intermediate solutions
        param_cg.xt=x;

        %% Preconditioning strategies

        [Ad,DA]=jacobi(Ar);

        for k=1:n_prec
            switch precond_tech{k}
                case 'No preconditioning'
                    [T(k).cond(i,j),lambda_min]=itcond(Ar);
                    [T(k).sol{i,j},T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Ar,b,tol*norm(b)*sqrt(lambda_min),maxiter,[],x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm(b);

                case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

                    [Ad,DA]=jacobi(Ar);
                    [T(k).cond(i,j),lambda_min]=itcond(Ad);
                    norm_b=normA(DA*DA,b);
                    [T(k).sol{i,j},T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Ar,b,tol*norm_b*sqrt(lambda_min),maxiter,@(x) DA*(DA*x),x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'SIPIC' % SIPIC preconditioner

                    [As,Ps,P2s,P1s]=sipic(Ar);
                    [T(k).cond(i,j),lambda_min]=itcond(As);
                    norm_b=normA(Ps*Ps',b);
                    [T(k).sol{i,j},T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Ar,b,tol*norm_b*sqrt(lambda_min),maxiter, @(x) Ps*(Ps'*x), x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'Schwarz' % Schwarz preconditioner

                    [Az,P1z,P2z,Pfun]=schwarz(sp_trimmed, msh_trimmed, Ar, int_dofs, T(k).param);
                    [T(k).cond(i,j),lambda_min]=itcond(Az);
                    norm_b=normA(@(x) DA*(Pfun(DA*x)),b);
                    [T(k).sol{i,j},T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = cg(Ar,b,tol*norm_b*sqrt(lambda_min),maxiter,@(x) DA*(Pfun(DA*x)),x0,param_cg);
                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;

                case 'Deflation' % Deflation based preconditioner

                    [Ap,P1p,P2p,Pfun,Ptfun,T(k).rank(i,j),~,ilM,Q]=deflation(sp_trimmed, msh_trimmed, Ar, int_dofs, T(k).param);
                    [T(k).cond(i,j),lambda_min]=itcond(Ap(ilM,ilM));
                    norm_b=normA(DA*DA,Pfun(b));
                    [T(k).sol{i,j},T(k).iter(i,j),T(k).res{i,j},T(k).res_prec{i,j},T(k).err{i,j}] = dpcg(Ar,b,tol*norm_b*sqrt(lambda_min),maxiter,Pfun,@(x) DA*(DA*x),x0,@(x) Ptfun(x)+Q(b),param_cg);

                    T(k).res_prec{i,j}=T(k).res_prec{i,j}/norm_b;


                otherwise
                    error('Not implemented')
            end
            T(k).err{i,j}=T(k).err{i,j}/normA(Ar,x);
            T(k).res{i,j}=T(k).res{i,j}/norm(b);
        end


    end

end

%% Post-processing
 
for i=1:ne
    for k=1:n_prec
        v=zeros(sp_cell{i}.ndof,T(k).iter(i,j)+1);
        v(int_dofs_cell{i},:)=T(k).sol{i,1};
        % Back-transform to the B-spline basis
        v=C_cell{i}*v;
        T(k).sol{i,1}=v;
        T(k).errl2{i,1}=sp_l2_error(sp_cell{i}, msh_cell{i}, T(k).sol{i,1}, uex);
    end
end

%% Condition number

% Condition number of the (preconditioned) mass
figure
for k=1:n_prec
    loglog(eta_vec, T(k).cond, 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
trendline(eta_vec(6:end), T(id_jac).cond(6:end), color=T(id_jac).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
trendline(eta_vec(10:20), T(id_sip).cond(10:20), color=T(id_sip).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
title('Conditioning - $K$')
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
title('Convergence history - $K$')
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
title('Convergence history - $K$')
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
title('Convergence history - $K$')
legend('Location', 'northeast')
legend show
xlabel('Iteration number')
ylabel('$A$-norm error')

% L2 Error
figure
for k=1:n_prec
    semilogy(T(k).errl2{s,1}, 'Linestyle', '-', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Convergence history - $K$')
legend('Location', 'northeast')
legend show
xlabel('Iteration number')
ylabel('$L^2$ Error')

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
legend(p, strcat({T.name}, {T.addendum}), 'Location', 'southwest')
legend show
ylim([1e-15 1e3])
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
legend(p, strcat({T.name}, {T.addendum}), 'Location', 'southwest')
legend show
ylim([1e-15 1e3])
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
legend(p, strcat({T.name}, {T.addendum}), 'Location', 'southwest')
legend show
ylim([1e-15 1e3])
xlabel('Iteration number')
ylabel('Relative error')

% L2 Error
figure
for i=1:ne
    for k=1:n_prec
        p(k)=semilogy(T(k).errl2{i,1}, 'Linestyle', '-', 'Marker', 'none', 'Color',  [T(k).color 0.5], 'LineWidth', 1);
        hold on; grid on;
    end
end
title('$L^2$ error')
legend(p, strcat({T.name}, {T.addendum}), 'Location', 'southwest')
legend show
ylim([1e-15 1e3])
xlabel('Iteration number')
ylabel('$L^2$ error')

%% Iteration numbers

figure
for k=1:n_prec
    semilogx(eta_vec, T(k).iter, 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Iteration counts - $K$')
legend('Location', 'northeast')
legend show
xlabel('$\eta$')
ylabel('Iteration count')

% Create table
tab=array2table([eta_vec [T.iter]], 'VariableNames', [{'eta'} precond_tech]);

%% Check eigenvalues

s=20;
K = op_gradu_gradv_trimming(sp_cell{s}, sp_cell{s}, msh_cell{s});
K=C_cell{s}'*K*C_cell{s};

Ar=K(int_dofs,int_dofs);
nk=500;

for k=1:n_prec

    switch precond_tech{k}
        case 'No preconditioning'

        case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

            Ad=jacobi(Ar);
            T(k).e=real(sort(eigs(Ad, nk, 'smallestabs')));

        case 'SIPIC' % SIPIC preconditioner

            As=sipic(Ar);
            T(k).e=real(sort(eigs(As, nk, 'smallestabs')));

        case 'Schwarz' % Schwarz preconditioner

            Az=schwarz(sp_cell{s}, msh_cell{s}, Ar, int_dofs, T(k).param);
            T(k).e=real(sort(eigs(Az, nk, 'smallestabs')));

        case 'Deflation' % Deflation based preconditioner

            [Ap,~,~,~,~,r,iss,isl]=deflation(sp_cell{s}, msh_cell{s}, Ar, int_dofs, T(k).param);


            rkr=nk-r;
            er=real(sort(eigs(Ap(isl,isl), rkr, 'smallestabs')));
            T(k).e=[zeros(r,1); er];

        otherwise
            error('Not implemented')

    end
end

figure
for k=1:n_prec
    semilogy(T(k).e(1:nk), 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
title('Spectrum of $K$')
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

for k=1:n_prec
    nexttile
    p(k)=semilogx(T(k).e, zeros(1,nk), 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum]);
    xlim([xmin xmax])
    grid on;
end

leg=legend(p, strcat({T.name}, {T.addendum}), 'Location', 'southoutside', 'NumColumns', n_prec);
leg.Layout.Tile = 'South';
sgtitle('Spectrum of $K$', 'Interpreter', 'Latex', 'FontSize', 14)
fontsize(fig, 14, 'points')

%% Rank values
ranks=cat(2,T(id_defl).rank);

s=1;
msh_plot_trim_geo(output{s})

[~,~,~,IS] = msh_split(sp_cell{s}, msh_cell{s}, 1);
ISr = linear_dep(sp_cell{s}, msh_cell{s}, 1);
msh_plot_supp(output{s}, sp_cell{s}, IS, 'red', true)
msh_plot_supp(output{s}, sp_cell{s}, ISr, 'blue', true)

%% L2 projection

rhs=op_f_v_gen(sp_cell{s}, msh_cell{s}, problem_data.uex);
M = op_u_v_trimming(sp_cell{s}, sp_cell{s}, msh_cell{s});
v_ref=M\rhs;

%% Visualize final solution
npts_per_el=[5,5];
[c, data_points, connectivity]=sp_eval(T(id_defl).sol{s,1}(:,end), sp_cell{s}, msh_cell{s}, geometry, npts_per_el);

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