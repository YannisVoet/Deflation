function [output, square_deg] = expanding_circle(eps_vec, d, ref)

% EXPANDING_CIRCLE: Creates a circular domain.
% INPUT:
% eps_vec:  circle radii
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
square = nrbsquare([0,0], lx, ly, 1, refinement);
square_deg = nrbsquare([0,0], lx, ly, d, refinement);

center = [0.5+1e-2, 0.5+1e-2];

for k=1:length(eps_vec)

    radius=eps_vec(k);

    % Trimming curve in the parametric domain
    % Circle
    loop_0 = struct();
    loop_0(1).curve = nrbcirc(radius, center);
    loop_0(1).label = 1;

    trimmed_srf.srf = square;
    trimmed_srf.trim_loops={loop_0};

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;

end

end