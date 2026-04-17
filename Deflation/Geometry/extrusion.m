function [output, square_deg] = extrusion(eps_vec, d, ref)

% EXTRUSION: Creates a plate with a cut-out.
% INPUT:
% eps_vec:  parameter controlling the radius of the cut-out
% d:        polynomial degree
% ref:      mesh refinement
% OUTPUT:
% output:   structure containing the trimmed surfaces
% square:   structure representing the fictitious domain

if nargin<3
    refinement=16; % for computations (the power of 2 should be greater
    % or equal to 3)
else
    refinement=ref;
end

% Construct cartesian grid to be cut
lx = 1;
ly = 1;
n_ref=0;
square = nrbsquare([0, 0], lx, ly, 1, refinement);
square_deg = nrbsquare([0,0], lx, ly, d, refinement);

for k = 1:length(eps_vec)

    eps = eps_vec(k);

    % Radius of the holes
    r=sqrt(5)/refinement(1)-eps;

    % Centers
    c5=[0.5,0.25];
    c7=[0.5,0.75];

    % Trimming reparametrization
    loop_1(1).curve = nrbline([0, 1], [1, 1]);
    loop_1(1).label = 4;
    loop_1(2).curve = nrbline([1, 1], [1, 0]);
    loop_1(2).label = 2;
    loop_1(3).curve = nrbline([1, 0], [0, 0]);
    loop_1(3).label = 3;
    loop_1(4).curve = nrbline([0, 0], [0, 1]);
    loop_1(4).label = 1;

    loop_2(1).curve = nrbcirc(r, c5, pi, 2*pi);
    loop_2(1).label = 5;
    loop_2(2).curve = nrbline([0.5+r, 0.25], [0.5+r, 0.75]);
    loop_2(2).label = 6;
    loop_2(3).curve = nrbcirc(r, c7, 0, pi);
    loop_2(3).label = 7;
    loop_2(4).curve = nrbline([0.5-r, 0.75], [0.5-r, 0.25]);
    loop_2(4).label = 8;

    all_loops = {loop_1, loop_2};

    trimmed_srf.srf=square;
    trimmed_srf.trim_loops=all_loops;

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;

end