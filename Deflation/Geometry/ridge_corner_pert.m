function [output, square_deg] = ridge_corner_pert(eps_vec, d, ref)

% CORNER_CUT_PERT: Perturbed "corner-cut" configuration.
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
corner=[0,0];

srf = nrb4surf (corner, corner+[lx 0], corner+[0 ly], corner+[lx ly]);

srf = nrbdegelev (srf, [0 0]);
srf_deg = nrbdegelev (srf, d-[1 1]);

[~,~,new_knots] = kntrefine (srf.knots{1}, refinement-1, 1, 0);

% Perturb the knots
seed=16; % Fix the random seed
rng(seed);
m=length(new_knots);
new_knotsp=new_knots+5e-2*rand(1,m);

square = nrbkntins (srf, {new_knotsp,new_knots});
square_deg = nrbkntins (srf_deg, {new_knotsp,new_knots});

for k = 1:length(eps_vec)

    eps = eps_vec(k);

    A = [0.1875, 0]; B = [0.1875, 0.5]; C = [0.5, 0.75+eps]; D = [0.8125, 0.5]; E = [0.8125, 0];
    loop_0 = struct();
    loop_0(1).curve = nrbline(A, B);
    loop_0(2).curve = nrbline(B, C);
    loop_0(3).curve = nrbline(C, D);
    loop_0(4).curve = nrbline(D, E);
    loop_0(5).curve = nrbline(E, A);
    loop_0(1).label = 5;
    loop_0(2).label = 6;
    loop_0(3).label = 7;
    loop_0(4).label = 8;
    loop_0(5).label = 9;
    trimmed_srf.srf=square;
    trimmed_srf.trim_loops={loop_0};

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;

end
