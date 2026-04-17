% Trimmed line [0, 0.75+eps]
% Plot the 1D spline and Lagrange bases, distinguishing between large, small 
% and inactive basis functions.

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

deltarange = [1e-2];

mrange = [3];
% mrange = [2];
drange = [2]; % degree

nat_frequencies = @(k, l_trimmed) (2*k - 1)/2 * pi/l_trimmed;
method_data.stabilization = false;
method_data.theta = 1;
method_data.stabilization_type = 'physical'; % 'parametric'

for d = drange
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
            method_data.degree     = d;           % Degree of the splines
            % method_data.regularity = d-1;
            method_data.regularity = 0;
            method_data.nquad      = d+1;

            method_data.solver      = 'jacobi';

            %% Manually creating the reparam structure
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
            for j = 1:numel(trim_elem_ids)
                trim_elems(j).elem_id = non_trim_elem_ids(numel(non_trim_elem_ids))+j;
                if( m == 1 )
                    % trim_elems(j).tiles = nrbsquare([0+(j-1)*element_size,0.5*length],element_size,0.25*length+delta,d,0);
                else
                    reparametrization = nrbline([0.75*lg,0],[0.75*lg + delta,0]);
                    trim_elems(j).tiles = nrbdegelev(reparametrization, d-1);
                    trim_elems(j).nb_tiles=1;
                    trim_elems(j).nb_pts=0;
                    trim_elems(j).quad_pts=[];
                    trim_elems(j).quad_weights=[];
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

            method_data.reparam = reparam;

            %% Solver
            % 3) CALL TO THE SOLVER=
            [geometry_grid, msh_trimmed, sp_trimmed, u, stiff_mat] = solve_laplace_trimming(problem_data, method_data);

            i = i+1;
        end
    end

    %% Plot spline basis

    [IL, IS] = msh_split(sp_trimmed, msh_trimmed);
    dofs=1:sp_trimmed.space_untrimmed.ndof;
    inact=setdiff(dofs, union(IL, IS));

    plotbasis(sp_trimmed.space_untrimmed, indices={IL, IS, inact}, linestyle={'-', ':', '-'}, linewidth={1, 2, 1}, colors={'b', 'b', [0.5 0.5 0.5]});
    xregion(reparam.domain(1), reparam.domain(2), FaceColor=[0.2 0.2 0.2]);
    ylim([0 1.5])
    title('Spline basis')
    box on

    %% Plot B-spline and Lagrange bases

    [~,C]=interpbasis(sp_trimmed.space_untrimmed, interp_pts='uniform', trunc=false);

    % Lagrange basis
    figure
    plotbasis(sp_trimmed.space_untrimmed,C=C);
    title('Lagrange basis')
    box on

    % B-spline basis
    figure
    plotbasis(sp_trimmed.space_untrimmed);
    title('B-spline basis')
    box on

end
