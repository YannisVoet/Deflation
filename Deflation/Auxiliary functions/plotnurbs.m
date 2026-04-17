function plotnurbs(nurbs, varargin)
%
% PLOTNURBS: plots a NURBS object (curve, surface, or the boundary of a
% NURBS volume).
%
% PLOTNURBS(nurbs) plots the NURBS object nurbs defined as a struct array.
% For multi-patch geometries, nurbs is a nonscalar struct array where
% nurbs(i) defines the ith patch. The geometry map for each patch is
% evaluated at 50 uniformily distributed points in each parametric direction.
% By default, black lines identify surface and volume edges.
%
% PLOTNURBS(nurbs, npts) specifies the number of points in each direction.
% If npts is a scalar, the same number of points is used along each
% parametric direction. Default: 50.
%
% PLOTNURBS(nurbs, npts, param) also provides optional parameters for
% changing the appearance. param is a struct array containing three fields:
%   'plot'     - plot properties for object boundaries (see properties of 
%                'plot3' or 'plot').
%   'surf'     - surface properties (see properties of 'surfl' or 'surf').
%   'axis'     - axis properties (see properties of 'axis').
% If nurbs is a nonscalar struct array (for a multi-patch geometry), param
% may also be a nonscalar struct array controlling the properties of each
% patch (the size of nurbs and param must then coincide). If param is a
% scalar struct array, the same properties apply to each patch. When 
% specifying optional parameters, leave npts empty '[]' to use default value.
%
% Some useful object properties are listed below:
% -'plot'
% 'LineWidth'
% 'Color'
%
% -'surf'
% 'FaceColor'
% 'EdgeColor'
% 'FaceAlpha'
%
% -'axis'-
% 'XLim'
% 'YLim'
% 'ZLim'
% 'CLim'
% 'XGrid' and 'XMinorGrid' (and the likes)
% 'View'
% 'Colormap'
% 'DataAspectRatio' Set to [1 1 1] for equal axes.
% 'XLabel.String'
% 'YLabel.String'
% 'ZLabel.String'
%
% For an exhaustive list, see MATLAB's documentation.
%
% Example 1: plot a multi-patch plate with a hole with 20 evaluation points
% along each direction for each patch.
%
% [geometry, ~,~,~,~] = mp_geo_load ('plate_with_hole_mp.txt');
% nurbs=[geometry.nurbs];
% plotnurbs(nurbs, 20);
% view(2)
%
% Example 2: plot a magnet shaped domain and set the linewidth of
% boundaries to 1.
%
% [geometry, ~,~,~,~] = mp_geo_load ('geo_magnet.txt');
% nurbs=[geometry.nurbs];
% param.plot.LineWidth=1;
% plotnurbs(nurbs, [], param);
%
% Example 3: plot a twisted 3-patch box and specify color and
% transparency of faces for each patch and linewidth of edges.
%
% [geometry, ~,~,~,~] = mp_geo_load ('geo_twisted_pipe_mp.txt');
% nurbs=[geometry.nurbs];
% 
% surf_prop=struct('FaceColor', {[0 0.3 0.7], [0 0.5 0.5], [0 0.7 0.3]}, 'FaceAlpha', 0.3);
% plot_prop=struct('LineWidth', 1.5);
% param=struct('plot', plot_prop, 'surf', {surf_prop(1), surf_prop(2), surf_prop(3)});
% 
% plotnurbs(nurbs, [], param);

%% Set default figure appearance

Sd.plot.Color='k';
Sd.surf.EdgeColor='none';

Sd.axis.Colormap=summer;
Sd.axis.DataAspectRatio=[1 1 1]; % 'equal' axes

%% Set algorithm parameters

% Set default parameters
Default{1}=50;
Default{2}=Sd;

% Replace empty inputs with default parameters
def=cell2mat(cellfun(@isempty, varargin, 'UniformOutput', false));
[varargin{def}]=Default{def};

Param = inputParser;
Param.addRequired('nurbs');
Param.addOptional('npts',   Default{1}, @(x) isvector(x) & all(x > 0));
Param.addOptional('param',  Default{2}, @isstruct);
Param.parse(nurbs, varargin{:});

%% Check algorithm parameters
npts=Param.Results.npts;
Se=Param.Results.param;


npatch=numel(nurbs);

if npatch>1 % Multi-patch geometry
    S=Sd;

    if numel(Se)==1
        for k=1:npatch
            S(k)=merge(Sd, Se);
        end
    else
        for k=1:npatch
            S(k)=merge(Sd, Se(k));
        end
    end

    for k=1:npatch
        plotnurbs(nurbs(k), npts, S(k));
        hold on
    end

    return

else
    S=merge(Sd, Se);
end


order=nurbs.order;
knots=nurbs.knots;

% Dimension
dim=length(order);

if length(npts)~=dim
    npts=npts*ones(1,dim);
end


switch dim
    case 1 % Plot a NURBS curve

        p = nrbeval(nurbs, linspace(knots(1), knots(end), npts));

        % 3D curve
        h=plot3(p(1,:), p(2,:), p(3,:));
        h=setObject(h, S.plot);
        grid on;


    case 2 % Plot a NURBS surface

        p = nrbeval(nurbs, {linspace(knots{1}(1),knots{1}(end),npts(1)) ...
            linspace(knots{2}(1),knots{2}(end),npts(2))});

        % light surface
        h=surfl(squeeze(p(1,:,:)), squeeze(p(2,:,:)), squeeze(p(3,:,:)));
        h=setObject(h, S.surf);
        hold on; grid on;


        % Plot the boundaries
        bnd = nrbextract(nurbs);

        plotnurbs(bnd(1), npts(2), S);
        plotnurbs(bnd(2), npts(2), S);
        plotnurbs(bnd(3), npts(1), S);
        plotnurbs(bnd(4), npts(1), S);


    case 3 % Plot the boundaries of a NURBS volume
        bnd = nrbextract(nurbs);

        plotnurbs(bnd(1), npts([2 3]), S);
        plotnurbs(bnd(2), npts([2 3]), S);
        plotnurbs(bnd(3), npts([1 3]), S);
        plotnurbs(bnd(4), npts([1 3]), S);
        plotnurbs(bnd(5), npts([1 2]), S);
        plotnurbs(bnd(6), npts([1 2]), S);
end

axis=gca;
axis=setObject(axis, S.axis);


%% Setting and updating graphic properties
    function[object]=setObject(object, param, varargin)

        if nargin>2
            X=varargin{1};
        else
            X={};
        end

        fields=setdiff(fieldnames(param), X, 'stable');
        s=numel(object);

        for i=1:length(fields)
            if isstruct(param.(fields{i}))
                object.(fields{i})=setObject(object.(fields{i}), param.(fields{i}));
            else

                if s==1
                    object.(fields{i})=param.(fields{i});
                else % Array of objects
                    if iscell(param.(fields{i})) % Objects take different parameters
                        [object.(fields{i})]=param.(fields{i}){:};
                    elseif size(param.(fields{i}),2)==1 % All objects take the same parameter
                        [object.(fields{i})]=deal(param.(fields{i}));
                    else % Objects take different parameters (stored in a regular array)
                        a=num2cell(param.(fields{i}),1);
                        [object.(fields{i})]=a{:};
                    end

                end
            end
        end

    end

end
