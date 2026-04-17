function [output, square_deg] = tips3(eps_vec, d, ref)

% TIPS3: Geometry with 3 tips for maximally smooth spline discretizations.
% INPUT:
% eps_vec:  parameter controlling the height of the ridge
% d:        polynomial degree
% ref:      mesh refinement
% OUTPUT:
% output:   structure containing the trimmed surfaces
% square:   structure representing the fictitious domain

% Construct cartesian grid to be cut
lx = 1;
ly = 1;

if nargin<3
    refinement=16;
else
    refinement=ref;
end

n_ref=0;
square = nrbsquare([0, 0], lx, ly, 1, refinement);
square_deg = nrbsquare([0,0], lx, ly, d, refinement);


for k = 1:length(eps_vec)

    eps = eps_vec(k);

    h=1/refinement(1);
    shift = 1/32;
    % shift = 1/25;
    A = [2*h, 0]; B = [2*h, 0.5]; C = [5.5*h, 12*h+eps]; D = [6.5*h, 11*h];
    E = [7.5*h 11*h]; F = [8.5*h 12*h+eps]; G = [9.5*h 11*h]; H = [10.5*h 11*h];
    I = [11.5*h 12*h+eps]; J = [15*h 0.5]; K = [15*h 0];

    loop_0 = struct();
    loop_0(1).curve = nrbline(A, B);
    loop_0(2).curve = nrbline(B, C);
    loop_0(3).curve = nrbline(C, D);
    loop_0(4).curve = nrbline(D, E);
    loop_0(5).curve = nrbline(E, F);
    loop_0(6).curve = nrbline(F, G);
    loop_0(7).curve = nrbline(G, H);
    loop_0(8).curve = nrbline(H, I);
    loop_0(9).curve = nrbline(I, J);
    loop_0(10).curve = nrbline(J, K);
    loop_0(11).curve = nrbline(K, A);
    loop_0(1).label = 5;
    loop_0(2).label = 6;
    loop_0(3).label = 7;
    loop_0(4).label = 8;
    loop_0(5).label = 9;
    loop_0(6).label = 10;
    loop_0(7).label = 11;
    loop_0(8).label = 12;
    loop_0(9).label = 13;
    loop_0(10).label = 14;
    loop_0(11).label = 15;
    trimmed_srf.srf=square;
    trimmed_srf.trim_loops={loop_0};

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;

end