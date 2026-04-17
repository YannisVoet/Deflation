function [output, square_deg] = square_with_hole(eps_vec, d, ref)

% SQUARE_WITH_HOLE: Rotated square with a hole in the middle (see [1]).
% INPUT:
% eps_vec:  rotation angles
% d:        polynomial degree
% ref:      mesh refinement
% OUTPUT:
% output:   structure containing the trimmed surfaces
% square:   structure representing the fictitious domain
%
% Reference: 
% [1] F. de Prenter, C. V. Verhoosel, G. J. van Zwieten, and E. H. van Brummelen. 
% Condition number analysis and preconditioning of the finite cell method. 
% Computer Methods in Applied Mechanics and Engineering, 316:297–327, 2017.

if nargin<3
    refinement=14;
else
    refinement=ref;
end

% Construct cartesian grid to be cut
lx = 1;
ly = 1;
n_ref=0;
square = nrbsquare([0, 0], lx, ly, 1, refinement);
square_deg = nrbsquare([0,0], lx, ly, d, refinement);

% square corners
s=0.3;
p1 = s*[1; -1];
p2 = s*[1; 1];
p3 = s*[-1; 1];
p4 = s*[-1; -1];

% geometry rotation
R = @(x) [cos(x) -sin(x); sin(x) cos(x)];
translation = 0.5;
rot=pi/4;

for k = 1:length(eps_vec)

    eps = eps_vec(k);

    % Radius of the holes
    r=sqrt(5)/refinement(1)-eps;

    % Centers
    c=[0.5,0.5];

    % Rotation and translation of square corners (rotations are done around the origin)
    p1_rot = R(rot) * p1 + translation;
    p2_rot = R(rot) * p2 + translation;
    p3_rot = R(rot) * p3 + translation;
    p4_rot = R(rot) * p4 + translation;

    % Trimming reparametrization
    loop_1(1).curve = nrbline(p1_rot, p2_rot);
    loop_1(1).label = 2;
    loop_1(2).curve = nrbline(p2_rot, p3_rot);
    loop_1(2).label = 4;
    loop_1(3).curve = nrbline(p3_rot, p4_rot);
    loop_1(3).label = 1;
    loop_1(4).curve = nrbline(p4_rot, p1_rot);
    loop_1(4).label = 3;

    loop_2(1).curve = nrbcirc(r, c, 0, 2*pi);
    loop_2(1).label = 5;

    all_loops = {loop_1, loop_2};

    trimmed_srf.srf=square;
    trimmed_srf.trim_loops=all_loops;

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;

end