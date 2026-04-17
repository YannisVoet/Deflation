% Verification of the scaling relations for the Jacobi preconditioner for 
% spline bases and various trimming configurations.

addpath(genpath('./Precond'))
addpath(genpath('./Auxiliary functions'))
addpath(genpath('./Geometry'))

clc
close all
clear variables

% Set the parameters
% polynomial degree
drange = [1:3];
% refinement level
n_ref = 0;
% case study
case_study = 'ridge'; % "middle-cut"
% case_study = 'ridge_corner'; % "corner-cut"
% case_study = 'translating_square_plate';
% case_study = 'ridge_corner_pert'; % perturbed "middle-cut"
eps_vec = logspace(-2,-4,20);
% eps_vec = logspace(-2,-3,20);

ne=length(eps_vec);
nd=length(drange);

% preconditioning techniques
precond_tech = {'Jacobi'};
n_prec=length(precond_tech);

% Load default parameters
T=def_param(precond_tech);

for k=1:n_prec
    T(k).cond_M=zeros(ne, nd);
    T(k).cond_K=zeros(ne, nd);
    T(k).rank=zeros(ne, nd);
end

% vector with the smallest volume ratio
eta_vec = zeros(ne,1);

for j = 1:length(drange)

    d = drange(j);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Degree %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', d);

    [output, square_deg] = feval(case_study,eps_vec,d);

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

        % Pure Neumann boundary conditions
        [K,M]=symmetrize(K,M);

        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp_trimmed, msh_trimmed, hfun, drchlt_sides);
        int_dofs = setdiff (1:sp_trimmed.ndof, drchlt_dofs);
        ndof=sp_trimmed.ndof;

        Kr = K(int_dofs, int_dofs);
        Mr = M(int_dofs, int_dofs);

        [Kr,Mr]=symmetrize(Kr,Mr);

        sK=1;
        sM=0;

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

                case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

                    [Kd,DK]=jacobi(Kr);
                    [Md,DM]=jacobi(Mr);

                    T(k).cond_K(i,j)=itcond(Kd,sK);
                    T(k).cond_M(i,j)=itcond(Md,sM);

                    [v1,eM1]=eigs(Mr, inv(DM^2), 1, 'smallestabs');

                    v1=v1/norm(v1);
                    T(k).r1(i,j)=v1'*Mr*v1;
                    T(k).r2(i,j)=v1'*((DM^2)\v1);

                case 'SIPIC' % SIPIC preconditioner
                    [Ks]=sipic(Kr);
                    [Ms]=sipic(Mr);

                    T(k).cond_K(i,j)=itcond(Ks,sK);
                    T(k).cond_M(i,j)=itcond(Ms,sM);

                case 'Schwarz' % Schwarz preconditioner

                    [Kz]=schwarz(sp_trimmed, msh_trimmed, Kr, int_dofs);
                    [Mz]=schwarz(sp_trimmed, msh_trimmed, Mr, int_dofs);

                    T(k).cond_K(i,j)=itcond(Kz,sK);
                    T(k).cond_M(i,j)=itcond(Mz,sM);

                case 'Deflation' % Deflation based preconditioner

                    [Kp,~,~,~,~,T(k).rank(i,j),~,ilK]=deflation(sp_trimmed, msh_trimmed, Kr, int_dofs, T(k).param);
                    [Mp,~,~,~,~,~,~,ilM]=deflation(sp_trimmed, msh_trimmed, Mr, int_dofs, T(k).param);

                    T(k).cond_K(i,j)=itcond(Kp(ilK,ilK),sK);
                    T(k).cond_M(i,j)=itcond(Mp(ilM,ilM),sM);

                otherwise
                    error('Not implemented')
            end
        end


    end

end

%% Plotting results

% Extract indices
id_def=cellfun(@(x) isequal(x,'No preconditioning'), precond_tech, 'UniformOutput', true);
id_jac=cellfun(@(x) isequal(x,'Jacobi'), precond_tech, 'UniformOutput', true);

c=colorgrad([0.5506 0.7515 0.8836], [0 0.1 0.9], nd);

% Condition number of the (preconditioned) mass
figure
for k=1:n_prec
    for j=1:nd
        loglog(eta_vec, T(k).cond_M(:,j), 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  c(j,:), 'LineWidth', 1, 'DisplayName', [precond_tech{k} ' ($p= ' num2str(drange(j)) '$)'])
        hold on; grid on;
    end

    for j=1:nd
        trendline(eta_vec, T(k).cond_M(:,j), color=c(j,:), linestyle='-', linewidth=0.5, variable='\eta', subset=true)
    end
end
title('Conditioning - $M$')
legend('Location', 'northeast', 'NumColumns', 2)
legend show
xlabel('$\eta$')
ylabel('$\kappa$')
ylim([1 1e20])

% Condition number of the (preconditioned) stiffness
figure
for k=1:n_prec
    for j=1:nd
        loglog(eta_vec, T(k).cond_K(:,j), 'Linestyle', 'none', 'Marker', T(k).marker, 'Color',  c(j,:), 'LineWidth', 1, 'DisplayName', [precond_tech{k} ' ($p= ' num2str(drange(j)) '$)'])
        hold on; grid on;
    end

    for j=1:nd
        trendline(eta_vec, T(k).cond_K(:,j), color=c(j,:), linestyle='-', linewidth=0.5, variable='\eta', subset=true)
    end
end
title('Conditioning - $K$')
legend('Location', 'northeast', 'NumColumns', 2)
legend show
xlabel('$\eta$')
ylabel('$\kappa$')
ylim([1 1e20])

%% Numerator and denominator of the Rayleigh quotient for Jacobi

figure
for j=1:nd
    loglog(eta_vec, T(id_jac).r1(:,j), 'Linestyle', 'none', 'Marker', 'v', 'Color',  c(j,:), 'LineWidth', 1, 'DisplayName', 'Numerator')
    hold on; grid on;
    loglog(eta_vec, T(id_jac).r2(:,j), 'Linestyle', 'none', 'Marker', T(id_jac).marker, 'Color',  c(j,:), 'LineWidth', 1, 'DisplayName', 'Denominator')
    trendline(eta_vec, T(id_jac).r1(:,j), color=c(j,:), linestyle='-', linewidth=0.5, variable='\eta', subset=true)
    trendline(eta_vec, T(id_jac).r2(:,j), color=c(j,:), linestyle='-', linewidth=0.5, variable='\eta', subset=true)
end
title('Scaling')
legend('Location', 'northeast')
legend show
xlabel('$\eta$')
ylabel('$\kappa$')
