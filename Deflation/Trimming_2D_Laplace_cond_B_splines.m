% Conditioning of the system matrices for the B-spline basis for various
% trimmed geometries and preconditioning strategies.

addpath(genpath('./Precond'))
addpath(genpath('./Auxiliary functions'))
addpath(genpath('./Geometry'))

clc
close all
clear variables

% Set the parameters
% refinement level
n_ref = 0;
% Geometry selection
geo=11;

switch geo
    case 1
        case_study = 'rotating_square';
        eps_vec = linspace(0, pi/4, 20);
        refinement = [16 16];
        drange = [3];
    case 2
        case_study = 'extrusion'; % For quadratic Lagrange bases
        eps_vec = logspace(-2,-4,20);
        refinement = [16 16];
        drange = [2];
    case 3
        case_study = 'tips2_C0';
        eps_vec = logspace(-2,-4,20);
        refinement = [16 16];
        drange = [3];
    case 4
        case_study = 'tips2';
        eps_vec = logspace(-2,-4,20);
        refinement = [16 16];
        drange = [3];
    case 5
        case_study = 'tips3';
        eps_vec = logspace(-2,-4,20);
        refinement = [16 16];
        drange = [3];
    case 6
        case_study = 'ridge';
        eps_vec = logspace(-2,-4,20);
        refinement = [16 16];
        drange = [3];
    case 7
        case_study = 'ridge_corner';
        eps_vec = logspace(-2,-4,20);
        refinement = [16 16];
        drange = [3];
    case 8
        case_study = 'square_with_hole';
        eps_vec = logspace(-2,-4,20);
        refinement = [16 16];
        drange = [3];
    case 9
        case_study = 'translating_square_plate';
        eps_vec = logspace(-2,-4,20);
        refinement = [16 16];
        drange = [3];
    case 10
        case_study = 'rotating_lattice_square';
        eps_vec = linspace(0, pi/4, 20);
        refinement = [16 16];
        drange = [3];
    case 11
        case_study = 'waveguide';
        eps_vec = logspace(-3,-4,20);
        refinement = [32 48];
        drange = [3];
end

ne=length(eps_vec);
nd=length(drange);

% preconditioning techniques
precond_tech = {'No preconditioning', 'Jacobi', 'SIPIC', 'Schwarz', 'Deflation'};
n_prec=length(precond_tech);

% Load default parameters
T=def_param(precond_tech);

for k=1:n_prec
    T(k).cond_M=zeros(ne, nd);
    T(k).cond_K=zeros(ne, nd);
    T(k).cond_A=zeros(ne, nd);
    T(k).rank=zeros(ne, nd);
end

id_def=cellfun(@(x) isequal(x,'No preconditioning'), {T.name}, 'UniformOutput', true);
id_jac=cellfun(@(x) isequal(x,'Jacobi'), {T.name}, 'UniformOutput', true);
id_sip=cellfun(@(x) isequal(x,'SIPIC'), {T.name}, 'UniformOutput', true);
id_sch=cellfun(@(x) isequal(x,'Schwarz'), {T.name}, 'UniformOutput', true);
id_defl=cellfun(@(x) isequal(x,'Deflation'), {T.name}, 'UniformOutput', true);

% vector with the smallest volume ratio
eta_vec = zeros(ne,1);

% Meshes and spaces
msh_cell=cell(ne,1);
sp_cell=cell(ne,1);

for j = 1:length(drange)

    d = drange(j);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Degree %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', d);

    [output, square_deg] = feval(case_study,eps_vec,d,refinement);

    problem_data.nmnn_sides = [];
    problem_data.drchlt_sides = [];
    problem_data.weak_drchlt_sides = [1 2 3 4];

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

        Kr = K(int_dofs, int_dofs);
        Mr = M(int_dofs, int_dofs);

        [Kr,Mr]=symmetrize(Kr,Mr);
        Ar=Kr+Mr;

        sK=1;
        sM=0;
        sA=0;

        %% Preconditioning strategies

        for k=1:n_prec
            switch precond_tech{k}
                case 'No preconditioning'

                    if i>2
                        c=polyfit(log(eta_vec(1:i-1)), log(eK1(1:i-1)), 1);
                    else
                        c=[1 1];
                    end

                    target=exp(c(2))*(eta_vec(i))^(2*degree(1));

                    if target < 1e-8
                        eigv=eigs(Kr, 10, target);
                        [~,idx]=min(abs(log(eigv)-log(target)));
                        eK1(i)=eigv(idx);
                    else
                        eK1(i)=max(eigs(Kr, 2, 'smallestabs'));
                    end

                    eKn(i)=max(eigs(Kr, 1, 'largestabs'));

                    T(k).cond_K(i,j)=eKn(i)/eK1(i);
                    T(k).cond_M(i,j)=itcond(Mr,sM);
                    T(k).cond_A(i,j)=itcond(Ar,sA);

                case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

                    [Kd,DK]=jacobi(Kr);
                    [Md,DM]=jacobi(Mr);
                    [Ad,DA]=jacobi(Ar);

                    T(k).cond_K(i,j)=itcond(Kd,sK);
                    T(k).cond_M(i,j)=itcond(Md,sM);
                    T(k).cond_A(i,j)=itcond(Ad,sA);

                    [v1,eM1]=eigs(Mr, inv(DM^2), 1, 'smallestabs');

                    v1=v1/norm(v1);
                    T(k).r1(i,j)=v1'*Mr*v1;
                    T(k).r2(i,j)=v1'*((DM^2)\v1);

                case 'SIPIC' % SIPIC preconditioner
                    [Ks]=sipic(Kr);
                    [Ms]=sipic(Mr);
                    [As]=sipic(Ar);

                    T(k).cond_K(i,j)=itcond(Ks,sK);
                    T(k).cond_M(i,j)=itcond(Ms,sM);
                    T(k).cond_A(i,j)=itcond(As,sA);

                case 'Schwarz' % Schwarz preconditioner

                    [Kz]=schwarz(sp_trimmed, msh_trimmed, Kr, int_dofs);
                    [Mz]=schwarz(sp_trimmed, msh_trimmed, Mr, int_dofs);
                    [Az]=schwarz(sp_trimmed, msh_trimmed, Ar, int_dofs);

                    T(k).cond_K(i,j)=itcond(Kz,sK);
                    T(k).cond_M(i,j)=itcond(Mz,sM);
                    T(k).cond_A(i,j)=itcond(Az,sA);

                case 'Deflation' % Deflation based preconditioner

                    [Kp,~,~,~,~,T(k).rank(i,j),~,ilK]=deflation(sp_trimmed, msh_trimmed, Kr, int_dofs, T(k).param);
                    [Mp,~,~,~,~,~,~,ilM]=deflation(sp_trimmed, msh_trimmed, Mr, int_dofs, T(k).param);
                    [Ap,~,~,~,~,~,~,ilA]=deflation(sp_trimmed, msh_trimmed, Ar, int_dofs, T(k).param);

                    T(k).cond_K(i,j)=itcond(Kp(ilK,ilK),sK);
                    T(k).cond_M(i,j)=itcond(Mp(ilM,ilM),sM);
                    T(k).cond_A(i,j)=itcond(Ap(ilA,ilA),sA);

                otherwise
                    error('Not implemented')
            end
        end


    end

end

%% Plotting results

% Condition number of the (preconditioned) mass
figure
for k=1:n_prec
    loglog(eta_vec, T(k).cond_M, 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
trendline(eta_vec, T(id_def).cond_M, color=T(id_def).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
trendline(eta_vec, T(id_jac).cond_M, color=T(id_jac).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
title('Conditioning - $M$')
legend('Location', 'northeast')
legend show
xlabel('$\eta$')
ylabel('$\kappa$')

% Condition number of the (preconditioned) stiffness
figure
for k=1:n_prec
    loglog(eta_vec, T(k).cond_K, 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
trendline(eta_vec, T(id_def).cond_K, color=T(id_def).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
trendline(eta_vec, T(id_jac).cond_K, color=T(id_jac).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
title('Conditioning - $K$')
legend('Location', 'northeast')
legend show
xlabel('$\eta$')
ylabel('$\kappa$')

% Condition number of the (preconditioned) linear combination stiffness + mass
figure
for k=1:n_prec
    loglog(eta_vec, T(k).cond_A, 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  T(k).color, 'LineWidth', 1, 'DisplayName', [precond_tech{k} T(k).addendum])
    hold on; grid on;
end
trendline(eta_vec, T(id_def).cond_A, color=T(id_def).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
trendline(eta_vec, T(id_jac).cond_A, color=T(id_jac).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
title('Conditioning - $K+M$')
legend('Location', 'northeast')
legend show
xlabel('$\eta$')
ylabel('$\kappa$')

%% Numerator and denominator of the Rayleigh quotient for Jacobi
figure
loglog(eta_vec, T(id_jac).r1, 'Linestyle', 'none', 'Marker', T(id_def).marker, 'Color',  T(id_def).color, 'LineWidth', 1, 'DisplayName', 'Numerator')
hold on; grid on;
loglog(eta_vec, T(id_jac).r2, 'Linestyle', 'none', 'Marker', T(id_jac).marker, 'Color',  T(id_jac).color, 'LineWidth', 1, 'DisplayName', 'Denominator')
trendline(eta_vec, T(id_jac).r1, color=T(id_def).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
trendline(eta_vec, T(id_jac).r2, color=T(id_jac).color, linestyle='-', linewidth=0.5, variable='\eta', subset=true)
title('Scaling')
legend('Location', 'northeast')
legend show
xlabel('$\eta$')
ylabel('$\kappa$')

%% Plot support of functions for the full and reduced set

s=1;
msh_plot_trim_geo(output{s})

[~,~,~,IS] = msh_split(sp_cell{s}, msh_cell{s}, 1);
ISr = linear_dep(sp_cell{s}, msh_cell{s}, 1);
msh_plot_supp(output{s}, sp_cell{s}, IS, 'red', true)
msh_plot_supp(output{s}, sp_cell{s}, ISr, 'blue', true)


