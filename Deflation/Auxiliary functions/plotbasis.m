function[]=plotbasis(space, options)

% PLOTBASIS: plots 1D spline basis functions (or their derivatives).
%
% PLOTBASIS(sp) plots the spline basis for the space sp.
%
% PLOTBASIS(sp,name,value) specifies optional name/value pair arguments:
%   'interval'  - Evaluation interval. Default: [0 1]
%   'r'         - Derivative order. Default: 0
%   'npts'      - Number of evaluation points. Default: 1000
%   'zpts'      - Abscissa. Leave empty ([]) to remove abscissa. Default: knots.
%   'C'         - Change of basis matrix. The new basis functions are
%                 defined as Bt = B*C, where Bt and B are row vectors.
%                 Default: identity.
%   'bfunc'     - Indices of basis functions to plot. Default: all indices.
%   'coeff'     - Plots the spline curve for the coefficient vector provided.
%   'indices'   - Cell array of groups of indices {G1,G2,...,Gm} for coloring.
%                 {G1,G2,...,Gm} must be pairwise disjoint. Indices that
%                 are in none of the groups Gi are treated as a default Matlab plot.
%                 Example: {1:space.ndof}. Default: empty.
%   'compact'   - Plots the basis functions on their support. This ensures 
%                 compactly supported functions are plotted as such.
%                 Default: true.
%   'colors'    - Cell array of colors {c1,c2,...,cm} (if indices are provided).
%                 Default: [0.5 0.5 0.5] (gray).
%   'linestyle' - Cell array of linestyles {l1,l2,...,lm} (if indices are provided).
%                 Default: '-' (continuous line).
%   'linewidth' - Cell array of linewidth {w1,w2,...,wm} (if indices are provided).
%                 Default: 0.5.
%
% Reference:
% [1] T. Lyche and K. Morken. Spline methods draft. Technical report, 
% Department of Mathematics, University of Oslo, 2018.

arguments
    space
    options.interval (1,2) {mustBeNumeric} = [0 1]
    options.r (1,1) {mustBeNumeric} = 0
    options.npts (1,1) {mustBeNumeric} = 1000
    options.zpts = space.knots{1}
    options.C = speye(space.ndof,space.ndof)
    options.bfunc = 1:space.ndof
    options.coeff = false
    options.indices = {}
    options.compact = true
    options.colors = {[0.5 0.5 0.5]}
    options.linestyle = {'-'}
    options.linewidth = {0.5}
end

% Knot vector
knots=space.knots{1};
% Degree
p=space.degree;
% Space dimension
n=space.ndof;
% Derivative
r=options.r;
% Evaluation interval
interval=options.interval;
% Number of points
npts=options.npts;
% Parametric points
x=linspace(interval(1), interval(2), npts);
% Change of basis matrix
C=options.C;

%% Ensure p+1 regularity
% This forms an augmented knot vector in case the knot vector is not 
% p+1 regular with p+1 common knots at the two ends (see Exercise 4.6 in [1])
d1=p+1-sum(knots==knots(1));
d2=p+1-sum(knots==knots(2));

knots_ext=sort([knots(1)*ones(1,d1) knots knots(end)*ones(1,d2)]);

% Dimension of the extended spline space (number of basis functions)
n_ext=length(knots_ext)-p-1;

% Knot span
mu=findspan_new(knots_ext,x,n_ext);
mu(mu>n_ext)=n;
% Evaluation of the basis functions
V=basisfunder_new(knots_ext, x, mu, p, r);

% Get indices
I=repmat((1:npts)', 1, p+1);
J=mu+(-p:0);

% The jth column contains the evaluation of the jth basis function
B=full(sparse(I(:), J(:), V(:), npts, n_ext));
% Remove artificial basis functions
B=B(:, d1+1:end-d2);
% Change basis
B=B*C;
options.bfunc = 1:size(B,2);
hold on

%% Plot basis

if isa(options.coeff, 'double')

    % Plot the spline curve
    plot(x, B*options.coeff, 'LineStyle', options.linestyle{1}, 'LineWidth', options.linewidth{1}, 'Color', options.colors{1})
else

    if options.compact
    B(B==0)=nan; % This ensures compactly supported functions are plotted as such
    end
    % Color the splines according to the group it belongs to
    for k=options.bfunc

        idx=cellfun(@(x) ismember(k,x), options.indices, 'UniformOutput', true);

        if sum(idx)>1
            error('An index cannot belong to multiple groups.')
        elseif sum(idx)==0
            plot(x, B(:,k)) % Default settings
        else
            plot(x, B(:,k), 'LineStyle', options.linestyle{idx}, 'LineWidth', options.linewidth{idx}, 'Color', options.colors{idx})
        end
    end
end

% Mark the abscissa
if ~isempty(options.zpts)
    z=zeros(1, length(options.zpts));
    plot(options.zpts, z, '.k', 'MarkerSize', 15)
end
end