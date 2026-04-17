% Trimmed line [0, 0.75+eps]
% Conditioning of the system matrices for the Lagrange basis and various
% preconditioning strategies.

addpath(genpath('./Precond'))
addpath(genpath('./Auxiliary functions'))
addpath(genpath('./Geometry'))

clear variables
close all
clc

% line geometry
lg = 1;
line = nrbline([0,0], [lg, 0]);
problem_data.geo_name = line;

% 1) PROBLEM DATA
problem_data.drchlt_sides = [1];
problem_data.nmnn_sides = [2];
problem_data.weak_drchlt_sides = [];

deltarange = logspace(-6, -2, 20);

mrange = [7]; % refinement level
drange = [3]; % degree

ne=length(deltarange);
nd=length(drange);

% preconditioning techniques
precond_tech = {'No preconditioning', 'Jacobi', 'SIPIC', 'Schwarz', 'Deflation'};
n_prec=length(precond_tech);

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

for j = 1:nd
    d = drange(j);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Degree %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',d);
    i = 1;

    for delta = deltarange
        n = 3;
        problem_data.f = @(x) n^2 * pi^2 * sin (n*pi*x);
        problem_data.hfun = @(x, ind) sin (n*pi*x);
        problem_data.g = @(x, ind) n * pi * cos (n*pi*x);
        problem_data.uex     = @(x) sin (n*pi*x);
        problem_data.graduex = @(x) n * pi * cos (n*pi*x);

        for m = mrange
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iteration %d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',i);

            h = 1 / (2.^m);

            method_data.nsub       = 2.^(m);
            method_data.degree     = d; % Degree of the splines
            method_data.regularity = 0;
            method_data.nquad      = d+1;

            method_data.solver      = 'jacobi';

            % Manually creating the reparam structure
            if( m == 1 ) % for now handling special case separately
                warning('You should start the simulation with at least 4 elements!')
                non_trimmed_elem_ids_cell{i} = [1];
                trimmed_ids_cell{i} = [2];
            else
                non_trimmed_elem_ids_cell{i} = [1:(2^m*3/4)];
                trimmed_ids_cell{i} = [2^m*3/4 + 1];
            end

            non_trim_elem_ids = non_trimmed_elem_ids_cell{i};
            trim_elem_ids = trimmed_ids_cell{i};

            nb_non_trim_elems=length(non_trim_elem_ids);
            nb_trim_elems=length(trim_elem_ids);

            % Creating the trim_elems structure
            for l = 1:numel(trim_elem_ids)
                trim_elems(l).elem_id = non_trim_elem_ids(numel(non_trim_elem_ids))+j;
                if( m == 1 )
                    % trim_elems(j).tiles = nrbsquare([0+(j-1)*element_size,0.5*length],element_size,0.25*length+delta,d,0);
                else
                    reparametrization = nrbline([0.75*lg,0],[0.75*lg + delta,0]);
                    trim_elems(l).tiles = nrbdegelev(reparametrization, d-1);
                    trim_elems(l).nb_tiles=1;
                    trim_elems(l).nb_pts=0;
                    trim_elems(l).quad_pts=[];
                    trim_elems(l).quad_weights=[];
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

            % Pure Neumann case
            int_dofs=1:ndof;

            Kr = K(int_dofs, int_dofs);
            Mr = M(int_dofs, int_dofs);

            [Kb,Mb]=symmetrize(Kr,Mr);

            % Expand matrices
            ndof_full=sp_trimmed.space_untrimmed.ndof;
            active_dofs=sp_trimmed.active_dofs;

            Mext=sparse(ndof_full, ndof_full);
            Kext=sparse(ndof_full, ndof_full);

            Mext(active_dofs, active_dofs)=Mb;
            Kext(active_dofs, active_dofs)=Kb;

            % Change basis from Bernstein to Lagrange polynomials on
            % uniform interpolation points.
            [Bext,Cext]=interpbasis(sp_trimmed.space_untrimmed, interp_pts='uniform', trunc=false);

            % System matrices with respect to the Lagrange basis
            Mext=Cext'*Mext*Cext;
            Kext=Cext'*Kext*Cext;

            % Restriction to active degrees of freedom
            M=Mext(active_dofs, active_dofs);
            K=Kext(active_dofs, active_dofs);

            Kr = K(int_dofs, int_dofs);
            Mr = M(int_dofs, int_dofs);

            [Kr,Mr]=symmetrize(Kr,Mr);

            C=Cext(active_dofs, active_dofs);
            B=Bext(active_dofs, active_dofs);

            sK=1;
            sM=0;

            %% Preconditioning strategies

            for k=1:n_prec
                switch precond_tech{k}
                    case 'No preconditioning'

                        if i>2
                            c=polyfit(log(deltarange(1:i-1)/h), log(eK1(1:i-1)), 1);
                        else
                            c=[1 1];
                        end

                        target=exp(c(2))*(deltarange(i)/h)^(2*degree-1);

                        if target < 1e-8
                            eigv=eigs(Kb, B'*B, 10, target);
                            [~,idx]=min(abs(log(eigv)-log(target)));
                            eK1(i)=eigv(idx);
                        else
                            eK1(i)=max(eigs(Kb, B'*B, 2, 'smallestabs'));
                        end

                        eKn(i)=max(eigs(Kb, B'*B, 1, 'largestabs'));

                        eM1=eigs(Mb, B'*B, 1, 'smallestabs');
                        eMn=eigs(Mb, B'*B, 1, 'largestabs');

                        T(k).cond_K(i,j)=eKn(i)/eK1(i);
                        T(k).cond_M(i,j)=eMn/eM1;

                    case 'Jacobi' % Diagonal scaling (Jacobi preconditioning)

                        [Kd,DK]=jacobi(Kr);
                        [Md,DM]=jacobi(Mr);


                        eKd1=max(eigs(Kb, B'*((DK^2)\B), 1+sK, 'smallestabs'));
                        eKdn=max(eigs(Kb, B'*((DK^2)\B), 1+sM, 'largestabs'));

                        eM1=eigs(Mb, B'*((DM^2)\B), 1, 'smallestabs');
                        eMn=eigs(Mb, B'*((DM^2)\B), 1, 'largestabs');

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

            tmp = msh_evaluate_element_list(msh_trimmed.msh_cart,[1]); % here we assume that all elements have the same size
            el_size(i) = max(tmp.element_size);
            dofs(i) = sp_trimmed.ndof;
            fprintf('Degree: %d. Total DOFs: %d \n', d, dofs(i));

            i = i+1;
        end
    end

    %% Plotting results

    % Plot basis
    figure
    plotbasis(sp_trimmed.space_untrimmed, C=Cext, zpts=[], npts=1e4)

    h=1/nsub;
    eta_vec = deltarange/h;

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

end
