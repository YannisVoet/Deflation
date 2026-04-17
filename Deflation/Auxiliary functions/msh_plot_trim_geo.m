function[]=msh_plot_trim_geo(output)

% MSH_PLOT_TRIM_GEO: Plots the trimmed geometry.

% Plot trimmed configuration
figure
srf = output.trim_srfs.srf;
nrbkntplot(srf, colormap=[0.9 0.9 0.9]); % Background grid
hold on; view(2)
trimmed_srfs_plot(output, 'color', [0.6 0.6 0.6]) % Physical domain
title('Fictitious and physical domain')

end