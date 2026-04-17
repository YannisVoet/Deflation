function [output, square_deg] = rotating_lattice_square(eps_vec, d, ref)

% ROTATING_LATTICE_SQUARE: Creates a lattice square with circular corner 
% exclusions and an interior circular exclusion (see [1]).
% INPUT:
% eps_vec:  rotation angles
% d:        polynomial degree
% ref:      mesh refinement
% OUTPUT:
% output:   structure containing the trimmed surfaces
% square:   structure representing the fictitious domain
%
% Reference:
% [1] F. de Prenter, C. Verhoosel, and E. Van Brummelen. Preconditioning 
% immersed isogeometric finite element methods with application to flow problems. 
% Computer Methods in Applied Mechanics and Engineering, 348:604–631, 2019.

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
% center of the interior circular exclusion
x_in = [9/20; 11/20];

% radii of the circular corner exclusions
R_ne = 1/3;
R_nw = 1/5;
R_sw = 1/pi;
R_se = 1/4;
% radius of the interior circular exclusion
R_in = 1/5;

% geometric transformation parameters
contr = 0.6; % contraction factor
transl = [0.2; 0.2];  % translation vector

% square corners after the geometric transformation
p1_map = contr*p1 + transl;
p2_map = contr*p2 + transl;
p3_map = contr*p3 + transl;
p4_map = contr*p4 + transl;
% center of the interior circular exclusion after the geometric transformation
x_in_map = contr*x_in + transl;

% radii of the circular corner exclusions after the geometric transformation
R_ne_map = contr*R_ne;
R_nw_map = contr*R_nw;
R_sw_map = contr*R_sw;
R_se_map = contr*R_se;
% radius of the interior circular exclusion after the geometric transformation
R_in_map = contr*R_in;

% auxiliary points
q1 = p1_map + [R_sw_map; 0];
q2 = p1_map + [0; R_sw_map];
q3 = p2_map - [0; R_nw_map];
q4 = p2_map + [R_nw_map; 0];
q5 = p3_map - [R_ne_map; 0];
q6 = p3_map - [0; R_ne_map];
q7 = p4_map + [0; R_se_map];
q8 = p4_map - [R_se_map; 0];

% geometry rotation
R = @(x) [cos(x) -sin(x); sin(x) cos(x)]; % rotation
% The center of the transformed geometry
origin = [(p1_map(1) + p4_map(1))/2 ; (p1_map(2) + p2_map(2))/2]; % center of the rotation

for k = 1:length(eps_vec)

    eps = eps_vec(k);

    % rotated square corners
    p1_rot = origin + R(eps)*(p1_map - origin);
    p2_rot = origin + R(eps)*(p2_map - origin);
    p3_rot = origin + R(eps)*(p3_map - origin);
    p4_rot = origin + R(eps)*(p4_map - origin);
    % rotated center of the interior circular exclusion
    x_in_rot = origin + R(eps)*(x_in_map - origin);

    % rotated auxiliary points
    q1_rot = origin + R(eps)*(q1 - origin);
    q2_rot = origin + R(eps)*(q2 - origin);
    q3_rot = origin + R(eps)*(q3 - origin);
    q4_rot = origin + R(eps)*(q4 - origin);
    q5_rot = origin + R(eps)*(q5 - origin);
    q6_rot = origin + R(eps)*(q6 - origin);
    q7_rot = origin + R(eps)*(q7 - origin);
    q8_rot = origin + R(eps)*(q8 - origin);

    % Trimming curve in the parametric domain
    loop_0 = struct();
    loop_0(1).curve = nrbcirc(R_sw_map, p1_rot, eps, pi/2+eps);
    loop_0(1).label = 1;
    loop_0(2).curve = nrbline(q2_rot, q3_rot);
    loop_0(2).label = 2;

    loop_0(3).curve = nrbcirc(R_nw_map, p2_rot, 3*pi/2+eps, 2*pi+eps);
    loop_0(3).label = 3;
    loop_0(4).curve = nrbline(q4_rot, q5_rot);
    loop_0(4).label = 4;

    loop_0(5).curve = nrbcirc(R_ne_map, p3_rot, pi+eps, 3*pi/2+eps);
    loop_0(5).label = 5;
    loop_0(6).curve = nrbline(q6_rot, q7_rot);
    loop_0(6).label = 6;

    loop_0(7).curve = nrbcirc(R_se_map, p4_rot, pi/2+eps, pi+eps);
    loop_0(7).label = 7;
    loop_0(8).curve = nrbline(q8_rot, q1_rot);
    loop_0(8).label = 8;

    loop_1(1).curve = nrbcirc(R_in_map, x_in_rot, eps, 2*pi+eps);
    loop_1(1).label = 9;

    trimmed_srf.srf=square;
    trimmed_srf.trim_loops={loop_0, loop_1};

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;
end
end