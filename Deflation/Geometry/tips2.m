function [output, square_deg] = tips2(eps_vec, d, ref)

% TIPS2: Geometry with 2 tips for maximally smooth spline discretizations.
% INPUT:
% eps_vec:  parameter controlling the height of the ridge
% d:        polynomial degree
% ref:      mesh refinement
% OUTPUT:
% output:   structure containing the trimmed surfaces
% square:   structure representing the fictitious domain

if nargin<3
    refinement=16;
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

    h=1/refinement(1);
    shift = 1/32;
    % shift = 1/25;
    A = [2*h+shift, 0]; B = [2*h+shift, 0.5]; C = [6*h+shift, 12*h+eps]; D = [7*h+shift, 11*h];
    E = [8*h+shift 11*h]; F = [9*h+shift 12*h+eps]; G = [14*h-shift 0.5]; H = [14*h-shift 0];

    loop_0 = struct();
    loop_0(1).curve = nrbline(A, B);
    loop_0(2).curve = nrbline(B, C);
    loop_0(3).curve = nrbline(C, D);
    loop_0(4).curve = nrbline(D, E);
    loop_0(5).curve = nrbline(E, F);
    loop_0(6).curve = nrbline(F, G);
    loop_0(7).curve = nrbline(G, H);
    loop_0(8).curve = nrbline(H, A);
    loop_0(1).label = 5;
    loop_0(2).label = 6;
    loop_0(3).label = 7;
    loop_0(4).label = 8;
    loop_0(5).label = 9;
    loop_0(6).label = 10;
    loop_0(7).label = 11;
    loop_0(8).label = 12;
    trimmed_srf.srf=square;
    trimmed_srf.trim_loops={loop_0};

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;

end