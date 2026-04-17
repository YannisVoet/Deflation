function[]=msh_plot_supp(output, sp_trimmed, I, color, plot_mesh)

% MSH_PLOT_SUPP: Plots the support of (trimmed) basis functions.
%
% MSH_PLOT_SUPP(output, sp_trimmed, I) plots the support of (trimmed) basis 
% functions for an index set I.
%
% MSH_PLOT_SUPP(output, sp_trimmed, I, color) specifies the color.
% Default: red.
%
% MSH_PLOT_SUPP(output, sp_trimmed, I, color, plot_mesh) specifies whether
% to plot the underlying mesh. Default: true.

arguments
    output
    sp_trimmed
    I
    color='red'
    plot_mesh=true
end

% Plot trimmed configuration
if plot_mesh
    figure
    srf = output.trim_srfs.srf;
    nrbkntplot(srf, colormap=[0.9 0.9 0.9]); % Background grid
    hold on;
    trimmed_srfs_plot(output, 'color', [0.6 0.6 0.6]) % Physical domain
end

m=length(I);
sp_untrimmed=sp_trimmed.space_untrimmed;
knots=sp_untrimmed.knots;
p=sp_untrimmed.degree;

for k=1:m
    % Index of global basis function
    ik=find(sp_trimmed.global_to_active==I(k));
    % Linear index
    [ix,iy]=ind2sub(sp_untrimmed.ndof_dir, ik);

    ix_in=knots{1}(ix); ix_end=knots{1}(ix+p(1)+1);
    iy_in=knots{2}(iy); iy_end=knots{2}(iy+p(2)+1);

    plot(ix_in, iy_in, '.', 'Color', color, 'MarkerSize', 15)
    supp(k)=nrbsquare([ix_in,iy_in] , ix_end-ix_in, iy_end-iy_in, 1);
end

if m>0
    param.surf.FaceColor='none';
    param.surf.EdgeColor='none';
    param.plot.Color=color;
    param.plot.LineWidth=1;
    param.axis.Colormap=[0.9 0.9 0.9];
    plotnurbs(supp, [], param) % Support
end

view(2)
title('Fictitious and physical domain')