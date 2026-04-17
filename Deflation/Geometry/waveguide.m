function [output, square_deg, c] = waveguide(delta_vec, d, ref)

% WAVEGUIDE: Spiky structure inspired from the waveguide geometry in 
% [1, Fig. 13]
% INPUT:
% eps_vec:  parameter controlling the height of the spikes
% d:        polynomial degree
% ref:      mesh refinement
% OUTPUT:
% output:   structure containing the trimmed surfaces
% square:   structure representing the fictitious domain
% c:        centers of the circular arcs
%
% Reference:
% [1] S. Eisenträger, L. Radtke, W. Garhuom, S. L¨ohnert, A. Düster, 
% D. Juhre, and D. Schillinger. An eigenvalue stabilization technique for 
% immersed boundary finite element methods in explicit dynamics.
% Computers & Mathematics with Applications, 166:129–168, 2024.

if nargin<3
    refinement=[16 16];
else
    if length(ref)==1
        refinement=[ref ref];
    else
        refinement=ref;
    end
end

% Construct cartesian grid to be cut
lx = 1;
ly = 1;
refx=refinement(1); hx=1/refx;
refy=refinement(2); hy=1/refy;
n_ref=0;
square = nrbsquare([0, 0], lx, ly, 1, refinement);
square_deg = nrbsquare([0,0], lx, ly, d, refinement);

nc=10;
c=zeros(nc,2);

for k = 1:length(delta_vec)

    shift=[hx/2 hy];

    delta = delta_vec(k);

    eta=delta;
    diam=d*hx-eta;
    r=0.5*diam;

    ltx=nc*(diam+eta);
    lty=0.5+delta;

    p1=[0 0]+shift;
    p2=[ltx 0]+shift;
    p3=[0 lty]+shift;
    p4=[ltx lty]+shift;

    loop_0 = struct();
    loop_0(1).curve = nrbline(p4, p2);
    loop_0(2).curve = nrbline(p2, p1);
    loop_0(3).curve = nrbline(p1, p3);
    loop_0(1).label = 2;
    loop_0(2).label = 3;
    loop_0(3).label = 1;

    p_old=p3+0.5*[eta 0];
    p_new=p_old+[diam 0];
    loop_0(4).curve = nrbline(p3, p_old);
    loop_0(4).label = 4;

    for j=1:nc-1
        c(j,:)=0.5*(p_old+p_new);
        loop_0(4+2*j-1).curve = nrbcirc(r, c(j,:), pi, 2*pi);
        loop_0(4+2*j-1).label=4+2*j-1;
        p_old=p_new+[eta 0];
        loop_0(4+2*j).curve = nrbline(p_new, p_old);
        loop_0(4+2*j).label=4+2*j;
        p_new=p_old+[diam 0];
    end

    c(nc,:)=0.5*(p_old+p_new);
    loop_0(4+2*nc-1).curve = nrbcirc(r, c(nc,:), pi, 2*pi);
    loop_0(4+2*nc-1).label=4+2*nc-1;
    p_old=p_new+0.5*[eta 0];
    loop_0(4+2*nc).curve = nrbline(p_new, p_old);
    loop_0(4+2*nc).label=4+2*nc;

    trimmed_srf.srf=square;
    trimmed_srf.trim_loops={loop_0};

    reparam = ref_trimmed_srfs(n_ref, trimmed_srf, 'reparam_deg', d);
    output{k} = reparam;

end