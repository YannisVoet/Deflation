function [output, square_deg] = ring(d, ref)

% RING: Creates a quarter of annulus (ring).
% INPUT:
% d:        polynomial degree
% ref:      mesh refinement
% OUTPUT:
% output:   structure containing the trimmed surfaces
% square:   structure representing the fictitious domain

%% This example crashes the trimmer

if nargin<3
    refinement=16; % for computations
else
    refinement=ref;
end

% Construct cartesian grid to be cut
lx = 1;
ly = 1;
n_ref=0;
square = nrbsquare([0, 0], lx, ly, 1, refinement);
square_deg = nrbsquare([0,0], lx, ly, d, refinement);

% Inner and outer radius
r1=0.5;
r2=1;

% Centers
c=[0 0];

% Trimming reparametrization
loop_1(1).curve = nrbline([0 r2], [0 r1]);
loop_1(1).label = 1;
loop_1(2).curve = nrbreverse(nrbcirc(r1, c, 0, pi/2));
loop_1(2).label = 3;
loop_1(3).curve = nrbline([r1, 0], [r2, 0]);
loop_1(3).label = 2;
loop_1(4).curve = nrbcirc(r2, c, 0, pi/2);
loop_1(4).label = 4;

all_loops = {loop_1};

trimmed_srf.srf=square;
trimmed_srf.trim_loops=all_loops;

reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
output{1} = reparam;

end