function [output, square_deg] = translating_rectangular_plate(eps_vec, d, ref)

% TRANSLATING_RECTANGULAR_PLATE: Trimmed rectangular square.
% INPUT:
% eps_vec:  parameter controlling the size of the smallest trimmed element
% d:        polynomial degree
% ref:      mesh refinement
% OUTPUT:
% output:   structure containing the trimmed surfaces
% square:   structure representing the fictitious domain

if nargin<3
    refinement=20;
else
    refinement=ref;
end

% Construct cartesian grid to be cut
lx = 1;
ly = 1;
n_ref=0;
square = nrbsquare([0, 0], lx, ly, 1, refinement);
square_deg = nrbsquare([0,0], lx, ly, d, refinement);

% square cornerns
p1 = [0; 0];
p2 = [0; 1];
p3 = [0.5; 1] ;
p4 = [0.5; 0];


for k = 1:length(eps_vec)

    eps = eps_vec(k);

    % Translation of 2 corners
    p3_tran = p3 + [eps; 0];
    p4_tran = p4 + [eps; 0];

    % Trimming curve in the parametric domain
    loop_0 = struct();
    loop_0(1).curve = nrbline(p1, p2);
    loop_0(2).curve = nrbline(p2, p3_tran);
    loop_0(3).curve = nrbline(p3_tran, p4_tran);
    loop_0(4).curve = nrbline(p4_tran, p1);

    loop_0(1).label = 1;
    loop_0(2).label = 2;
    loop_0(3).label = 3;
    loop_0(4).label = 4;
    trimmed_srf.srf=square;
    trimmed_srf.trim_loops={loop_0};

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;

end