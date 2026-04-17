% Experiment 3: Wave equation on a waveguide-inspired spiky structure.
% Complete simulation.

addpath(genpath('./Precond'))
addpath(genpath('./Auxiliary functions'))
addpath(genpath('./Geometry'))

clc
close all
clear variables

% Set the parameters
% Time span
I=[0 0.4];
% polynomial degree
drange = [3];
% refinement level
refinement = [32 48];
n_ref = 0;
% case study
case_study = 'waveguide';
eps_vec = logspace(-2,-3,20);

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


j=1;
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

i = 5;
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

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh_trimmed = msh_trimming (zeta, qn, qw, geometry, reparam);

% Construct space structure
sp_trimmed = sp_trimming (knots, degree, msh_trimmed);

% Assemble the matrices
K = op_gradu_gradv_trimming (sp_trimmed, sp_trimmed, msh_trimmed);
M = op_u_v_trimming (sp_trimmed, sp_trimmed, msh_trimmed);

[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp_trimmed, msh_trimmed, @(x, y, ind) zeros(size(x)), drchlt_sides);
int_dofs = setdiff (1:sp_trimmed.ndof, drchlt_dofs);
n=length(int_dofs);

Kr=K(int_dofs, int_dofs);
Mr=M(int_dofs, int_dofs);

[Kr,Mr]=symmetrize(Kr,Mr);

Nf=int_dofs;
Nd=drchlt_dofs;

nf=length(Nf);
nd=length(Nd);
nb_dof=sp_trimmed.ndof;

Data.Nf=Nf;
Data.Nd=Nd;

Data.sp=sp_trimmed;
Data.msh=msh_trimmed;
Data.geometry=geometry;
Data.nmnn_sides=nmnn_sides;
Data.drchlt_sides=drchlt_sides;

n_sp=sp_trimmed.space_untrimmed.ndof;
active_dofs=sp_trimmed.active_dofs;

% Number of sub-intervals in time
N=600;

% Implicit unconditionally stable method
Options=struct('param', [1/4 1/2], 'I', I, 'N', N, 'V', spdiags(1./sqrt(diag(Mr+Kr)), 0, n, n));

%% Computation of the numerical solution

[Data, Solution]=solveode([], Data, M, sparse(nb_dof, nb_dof), K, Options);

%% Post-processing

U=Solution.U;
d1_U=Solution.d1_U;

Discretization=[Data.Discretization];

N=Discretization.N;
T=Discretization.time_vec;

npts_per_el=[5,5];

[c, data_points, connectivity]=sp_eval(U, sp_trimmed, msh_trimmed, geometry, npts_per_el);

xp=data_points(:,1);
yp=data_points(:,2);

%% Numerical solution
Max=max(c, [], 'all');
Min=min(c, [], 'all');

clear S F
% Figure appearance
F.FigParam.Units='normalized';
F.FigParam.Position=[0 0 0.3 0.3];

S.type='patch';
S.PlotParam=struct('Vertices', data_points, 'Faces', connectivity);
S.AxisParam.CLim=[Min Max];
S.AxisParam.DataAspectRatio=[1 1 1];
S.AxisParam.Title.String='Solution';
S.DynamicParam=struct('CData', num2cell(c, 1));
S.ColorbarParam.Visible='on';
S.step=1;

F.view_video=false;
F.write_video=false;
F.shrink=true;
F.vid_name='Trimming_2D_Laplace_rotating_square';

plotsol(F,S)


%% Snapshots

% Number of snapshots
ns=4;

clear S F
% Figure appearance
F.FigParam.Units='normalized';
F.FigParam.Position=[0 0 1 0.4];
F.TiledParam.GridSize=[1 ns];
F.TiledParam.TileSpacing='compact';
F.TiledParam.Padding='none';

S.type='patch';
S.PlotParam=struct('Vertices', data_points, 'Faces', connectivity);
S.AxisParam.CLim=[Min Max];
% S.AxisParam.DataAspectRatio=[1 1 1];
S.AxisParam.FontSize=22;
S.AxisParam.Title.String='Solution';
S.DynamicParam=struct('CData', num2cell(c, 1));
S.ColorbarParam.Visible='off';
S.step=1;

% Copying
S(1:ns)=S;

steps=num2cell(round(linspace(1,N+1,ns)));
[S.step]=steps{:};

S(end).ColorbarParam.Visible='on';

F.view_video=false;
F.write_video=false;
F.shrink=true;
F.vid_name='1D_Laplace_trimming_exact_solution';

plotsol(F,S)
