function [output, square_deg] = rotating_square_circular_hole(eps_vec, d, ref)

% ROTATING_SQUARE_CIRCULAR_HOLE: Creates a rotated square with a circular 
% corner exclusion. The rotation is around this corner (see [1]).
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
p1 = [0; 0];
p2 = [0; 1];
p3 = [1; 1] ;
p4 = [1; 0];

% radii of the circular corner exclusions
R_sw = 3/(2*pi);

% geometric transformation parameters
contr = 0.31; % contraction factor
transl = [0.51; 0.51];  % translation vector

% square corners after the geometric transformation
p1_map = contr*p1 + transl;
p2_map = contr*p2 + transl;
p3_map = contr*p3 + transl;
p4_map = contr*p4 + transl;

% radii of the circular corner exclusion after the geometric transformation
R_sw_map = contr*R_sw;

% auxiliary points
q1 = p1_map + [R_sw_map; 0];
q2 = p1_map + [0; R_sw_map];

% geometry rotation
R = @(x) [cos(x) -sin(x); sin(x) cos(x)]; % rotation
% The center of the transformed geometry
origin = p1_map; % center of the rotation

for k = 1:length(eps_vec)

    eps = eps_vec(k);

    % rotated square corners
    p1_rot = origin + R(eps)*(p1_map - origin);
    p2_rot = origin + R(eps)*(p2_map - origin);
    p3_rot = origin + R(eps)*(p3_map - origin);
    p4_rot = origin + R(eps)*(p4_map - origin);

    % rotated auxiliary points
    q1_rot = origin + R(eps)*(q1 - origin);
    q2_rot = origin + R(eps)*(q2 - origin);

    % Trimming curve in the parametric domain
    loop_0 = struct();
    loop_0(1).curve = nrbcirc(R_sw_map, p1_rot, eps, pi/2+eps);
    loop_0(1).label = 1;

    loop_0(2).curve = nrbline(q2_rot, p2_rot);
    loop_0(2).label = 2;

    loop_0(3).curve = nrbline(p2_rot, p3_rot);
    loop_0(3).label = 3;

    loop_0(4).curve = nrbline(p3_rot, p4_rot);
    loop_0(4).label = 4;

    loop_0(5).curve = nrbline(p4_rot, q1_rot);
    loop_0(5).label = 5;

    trimmed_srf.srf=square;
    trimmed_srf.trim_loops={loop_0};

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;
end
end