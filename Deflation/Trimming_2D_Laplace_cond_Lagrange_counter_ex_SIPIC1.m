% First counter-example for SIPIC with the Lagrange basis.

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
refinement = 20;
n_ref = 0;
case_study = 'extrusion';
eps_vec = logspace(-2,-4,20);

ne=length(eps_vec);
nd=length(drange);

% preconditioning techniques
precond_tech = {'No preconditioning', 'Jacobi', 'SIPIC', 'Schwarz', 'Deflation'};
n_prec=length(precond_tech);
% vector with the smallest volume ratio
eta_vec = zeros(ne,1);

% Load default parameters
T=def_param(precond_tech);

for k=1:n_prec
    T(k).cond_M=zeros(ne, nd);
    T(k).cond_K=zeros(ne, nd);
    T(k).rank=zeros(ne, nd);
end

id_def=cellfun(@(x) isequal(x,'No preconditioning'), {T.name}, 'UniformOutput', true);
id_jac=cellfun(@(x) isequal(x,'Jacobi'), {T.name}, 'UniformOutput', true);
id_sip=cellfun(@(x) isequal(x,'SIPIC'), {T.name}, 'UniformOutput', true);
id_sch=cellfun(@(x) isequal(x,'Schwarz'), {T.name}, 'UniformOutput', true);
id_defl=cellfun(@(x) isequal(x,'Deflation'), {T.name}, 'UniformOutput', true);

for j = 1:length(drange)

    d = drange(j);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Degree %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n', d);

    [output, square_deg] = feval(case_study,eps_vec,d);

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
    method_data.regularity = [0 0];     % Regularity of the splines
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
        [Kb,Mb]=symmetrize(K,M);

        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp_trimmed, msh_trimmed, hfun, drchlt_sides);
        int_dofs = setdiff (1:sp_trimmed.ndof, drchlt_dofs);
        ndof=sp_trimmed.ndof;

        % Expand matrices
        ndof_full=sp_trimmed.space_untrimmed.ndof;
        active_dofs=sp_trimmed.active_dofs;

        Mext=sparse(ndof_full, ndof_full);
        Kext=sparse(ndof_full, ndof_full);

        Mext(active_dofs, active_dofs)=Mb;
        Kext(active_dofs, active_dofs)=Kb;

        % Change basis from Bernstein to Lagrange polynomials on
        % uniform interpolation points.
        [B,C]=interpbasis(sp_trimmed.space_untrimmed, interp_pts='uniform', trunc=false);

        % System matrices with respect to the Lagrange basis
        Mext=C'*Mext*C;
        Kext=C'*Kext*C;

        % Restriction to active degrees of freedom
        M=Mext(active_dofs, active_dofs);
        K=Kext(active_dofs, active_dofs);

        Kr = K(int_dofs, int_dofs);
        Mr = M(int_dofs, int_dofs);

        [Kr,Mr]=symmetrize(Kr,Mr);

        C=C(active_dofs, active_dofs);
        B=B(active_dofs, active_dofs);

        sK=1;
        sM=0;

        [~,DKb]=jacobi(Kb);
        [~,DMb]=jacobi(Mb);

        %% Preconditioning strategies

        for k=1:n_prec
            switch precond_tech{k}
                case 'No preconditioning'

                    if i>2
                        c=polyfit(log(eta_vec(1:i-1)), log(eK1(1:i-1)), 1);
                    else
                        c=[1 1];
                    end

                    target=exp(c(2))*eta_vec(i)^(2*degree(1));

                    if target < 1e-8
                        eigv=eigs(Kb, B'*B, 10, target);
                        [~,idx]=min(abs(log(eigv)-log(target)));
                        eK1(i)=eigv(idx);
                    else
                        eK1(i)=max(eigs(DKb*Kb*DKb, DKb*(B'*B)*DKb, 2, 'smallestabs'));
                    end

                    eKn(i)=max(eigs(Kb, B'*B, 1, 'largestabs'));

                    eM1=eigs(DMb*Mb*DMb, DMb*(B'*B)*DMb, 1, 'smallestabs');
                    eMn=eigs(Mb, B'*B, 1, 'largestabs');

                    T(k).cond_K(i,j)=eKn(i)/eK1(i);
                    T(k).cond_M(i,j)=eMn/eM1;

                case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

                    [Kd,DK]=jacobi(Kr);
                    [Md,DM]=jacobi(Mr);


                    eKd1=max(eigs(DKb*Kb*DKb, DKb*B'*((DK^2)\B)*DKb, 1+sK, 'smallestabs'));
                    eKdn=max(eigs(Kb, B'*((DK^2)\B), 1+sM, 'largestabs'));

                    [v1,eM1]=eigs(Mb, B'*((DM^2)\B), 1, 'smallestabs');
                    eMn=eigs(Mb, B'*((DM^2)\B), 1, 'largestabs');

                    v1=v1/norm(v1);
                    T(k).r1(i,j)=v1'*Mb*v1;
                    T(k).r2(i,j)=v1'*B'*((DM^2)\B)*v1;


                    T(k).cond_K(i,j)=eKdn/eKd1;
                    T(k).cond_M(i,j)=eMn/eM1;

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

