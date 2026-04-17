function [output, square_deg] = rotating_square(eps_vec, d, ref)

% ROTATING_SQUARE: Creates a rotated square (without any holes or cut-outs)
% INPUT:
% eps_vec:  rotation angles
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

% square corners
s=0.2+1e-5;
p1 = s*[1; -1];
p2 = s*[1; 1];
p3 = s*[-1; 1];
p4 = s*[-1; -1];

% geometry rotation
R = @(x) [cos(x) -sin(x); sin(x) cos(x)];
translation = 0.51;


for k = 1:length(eps_vec)

    eps = eps_vec(k);

    % Rotation and translation of square corners (rotations are done around the origin)
    p1_rot = R(eps) * p1 + translation;
    p2_rot = R(eps) * p2 + translation;
    p3_rot = R(eps) * p3 + translation;
    p4_rot = R(eps) * p4 + translation;

    % Trimming curve in the parametric domain
    loop_0 = struct();
    loop_0(1).curve = nrbline(p1_rot, p2_rot);
    loop_0(2).curve = nrbline(p2_rot, p3_rot);
    loop_0(3).curve = nrbline(p3_rot, p4_rot);
    loop_0(4).curve = nrbline(p4_rot, p1_rot);

    loop_0(1).label = 1;
    loop_0(2).label = 2;
    loop_0(3).label = 3;
    loop_0(4).label = 4;
    trimmed_srf.srf=square;
    trimmed_srf.trim_loops={loop_0};

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;

end