% Rotated square example: visualization of the active mesh and partitioning 
% into cut and uncut regions.

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
n_ref = 0;
% case study
case_study = 'rotating_square';
eps_vec = linspace(0, pi/4, 20);

ne=length(eps_vec);
nd=length(drange);

% Meshes and spaces
msh_cell=cell(ne,1);
sp_cell=cell(ne,1);

% vector with the smallest volume ratio
eta_vec = zeros(ne,1);

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

        % Save space and mesh
        msh_cell{i}=msh_trimmed;
        sp_cell{i}=sp_trimmed;

    end

end

%% Plot cut and uncut elements

s=7;

% Plot active mesh
msh_plot_trim_geo(output{s})
[~,~,~,~, good_ids, bad_trimmed_ids] = msh_split(sp_cell{s}, msh_cell{s}, 1);
msh_plot_elem(output{s}, msh_cell{s}, union(good_ids, bad_trimmed_ids), [0.9 0.5 0.1])
title('Active mesh')

% Plot cut and uncut elements
msh_plot_trim_geo(output{s})
msh_plot_elem(output{s}, msh_cell{s}, {good_ids, bad_trimmed_ids}, {[0 0.2 0.8], [0.9 0.1 0.1]})
title('Mesh partitioning')

