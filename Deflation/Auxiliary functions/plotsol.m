function[]=plotsol(F, varargin)

% PLOTSOL: plots the solution.
%
% PLOTSOL(F,S1,S2,...,Sn) creates a figure with n tiles and plots in the
% kth tile the contents of the kth structure Sk for k = 1,...,n.
%
% F is a structure containing the following fields:
%   'view_video'    - boolean parameter for viewing the video of a
%                     dynamical simulation. Default: false.
%   'write_video'   - boolean parameter for saving the video to MP4 format.
%                     Default: false.
%   'vid_duration'  - video duration (in seconds). Default: 10 s.
%   'vid_name'      - video name. Default: 'Output.mp4'.
%   'shrink'        - controls how simulations for different numbers of
%                     time steps are matched together. If 'shrink' is
%                     true, the number of visualization steps is equal to
%                     the smallest number of time steps. If 'shrink' is
%                     false, the number of visualization steps is equal to
%                     the largest number of time steps. Default: true.
%   'FigParam'      - figure properties (see properties of 'figure').
%   'TiledParam'    - tiledlayout properties (see properties of
%                     'tiledlayout').
%   'LegendParam'   - figure legend properties (see properties of 'legend').
%
% The (nonscalar) structures Sk contain (at most) the following fields:
%   'type'          - plot type (e.g 'patch', 'surf', etc).
%   'step_count'    - Displays the time step countdown (for
%                     dynamical simulations). Default: true.
%   'step'          - time step at which the solution is viewed (for
%                     dynamical simulations). Default: 1.
%   'PlotParam'     - plot properties (see properties of 'patch', 'surf',
%                     'quiver', etc).
%   'AxisParam'     - axis properties (see properties of 'axis').
%   'LegendParam'   - legend properties (see properties of 'legend').
%   'ColorbarParam' - colorbar properties (see properties of 'colorbar').
%   'DynamicParam'  - properties of 'PlotParam' that change during the
%                     course of a dynamical simulation.
%
% Some useful object properties are listed below:
% -'figure'
% 'Units'
% 'Position'
%
% -'tiledlayout'
% 'TileSpacing'
% 'Padding'
% 'GridSize'
%
% -'patch'-
% 'Vertices'
% 'Faces'
% 'CData' or 'FaceVertexCData'
% 'FaceColor'
% 'EdgeColor'
% 'LineWidth'
%
% -'quiver'-
% 'XData'
% 'YData'
% 'UData'
% 'VData'
% 'LineWidth'
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
% -'legend'-
% 'Visible'
% 'Location'
% 'Title.String'
% 'Layout.Tile'
%
% -'colorbar'-
% 'Visible'
% 'Label.String'
% 'Label.FontSize'
%
% For an exhaustive list, see MATLAB's documentation.
%
% Example 1: View the solution at time step 15000 of a dynamical finite
% element simulation with a 'patch' plot. The simulation data is stored
% along the columns of a matrix c.
%
% F.view_video=false;
%
% S.type='patch';
% S.PlotParam=struct('Vertices', vertices, 'Faces', faces);
% S.AxisParam.Title.String='Solution';
% S.DynamicParam=struct('CData', num2cell(c, 1));
% S.step=15000;
%
% plotsol(F,S)
%
% Example 2: View the simulations for input data c1 and c2 over the same
% mesh.
%
% F.view_video=true;
% F.shink=true;
%
% S1.type='patch';
% S1.PlotParam=struct('Vertices', vertices, 'Faces', faces);
% S1.AxisParam.Title.String='Solution 1';
% S1.DynamicParam=struct('CData', num2cell(c1, 1));
%
% S2=S1;
% S2.AxisParam.Title.String='Solution 2';
% S2.DynamicParam=struct('CData', num2cell(c2, 1));
%
% plotsol(F,S1,S2)
%
% Example 3: View a finite element mesh.
%
% F=struct([]);
%
% S.type='patch';
% S.PlotParam=struct('Vertices', vertices, 'Faces', faces);
% S.PlotParam.EdgeColor='black';
% S.PlotParam.FaceColor='none';
% S.AxisParam.Title.String='Mesh';
%
% plotsol(F,S)
%
% Example 4: Plot the results of three 1D dynamical simulations in the
% same figure. The data is stored along the columns of arrays c1, c2 and c3.
%
% c=cat(3, c1, c2, c3);
% c=permute(c, [1 3 2]);
% c=num2cell(c, [1 2]);
% c=squeeze(c);
%
% F.view_video=true;
% F.write_video=true;
% F.shrink=true;
% F.vid_name='1D_Dynamics';
%
% S.type='plot';
% S.PlotParam.XData=x;
% S.PlotParam.Color={'r', 'k', 'b'};
% S.AxisParam.YLim=[Min Max];
% S.AxisParam.Title.String='Solutions';
% S.LegendParam.Visible='on';
% S.LegendParam.Title.String={'Solution 1', 'Solution 2', 'Solution 3'};
% S.DynamicParam=struct('YData', c);
%
% plotsol(F,S)
%
% Example 5: View the results of a 2D multipatch dynamical simulation. eu
% is a cell array with as many elements as patches. Each element stores an
% np x np x (N + 1) tensor, where np is the number of evaluation points in
% each direction and N is the number of subdivisions in time.
%
% eu=cat(4,eu{:});
% eu=squeeze(num2cell(eu, [1 2 4]));
% eu=cellfun(@(x) squeeze(num2cell(squeeze(x), [1 2])), eu, 'UniformOutput', false);
%
% F.view_video=true;
% F.write_video=true;
% F.vid_name='2D_Multipatch';
%
% S.type='surf';
% S.PlotParam.XData=X;
% S.PlotParam.YData=Y;
% S.DynamicParam=struct('ZData', eu);
% S.AxisParam.DataAspectRatio=[1 1 1];
% S.AxisParam.CLim=[Min Max];
% S.AxisParam.Title.String='Solution';
%
% plotsol(F,S);
%
% Example 6: Plot solution snapshots at 6 different times with a single
% legend for the entire figure. This example illustrates the use of
% nonscalar struct arrays.
%
% F.view_video=false;
% F.write_video=false;
% F.shrink=true;
% F.vid_name='1D_Snapshots';
%
% F.FigParam.Units='normalized';
% F.FigParam.Position=[0 0 0.6 0.4];
% F.TiledParam.GridSize=[2 3];
% F.LegendParam.Visible='on';
% F.LegendParam.Layout.Tile='south';
% F.LegendParam.Title.String={'Solution 1', 'Solution 2', 'Solution 3'};
%
% S.type='plot';
% S.PlotParam.XData=x;
% S.PlotParam.Color={'r', 'k', 'b'};
% S.AxisParam.YLim=[Min Max];
% S.AxisParam.Title.String='Solutions';
% S.DynamicParam=struct('YData', c);
%
% S(1:6)=S;
%
% steps=num2cell(round(1/6*(1:6)*(N+1)));
% [S.step]=steps{:};
%
% plotsol(F,S)

%% Set default fields

% Set default video parameters
Fd.view_video=false;
Fd.write_video=false;
Fd.vid_duration=10; % seconds
Fd.vid_name='Output.mp4';
Fd.shrink=true;

% Set default figure parameters
Fd.FigParam.Units='normalized';
Fd.FigParam.Position=[0 1 1/2 1/2];

% Set default tiledlayout parameters
Fd.TiledParam.TileSpacing='tight';
Fd.TiledParam.Padding='tight';

% Set default figure legend parameters
Fd.LegendParam.Title.String={};
Fd.LegendParam.Visible='off';


struct_fields={'type', 'step_count', 'step', 'PlotParam', 'AxisParam', ...
    'LegendParam', 'ColorbarParam', 'DynamicParam'};

% Set default plot parameters
% Standard 2D plot
Sd.PlotParam.plot.LineWidth=1.5;

% Standard loglog plot
Sd.PlotParam.loglog.LineWidth=1.5;

% Standard semilogx plot
Sd.PlotParam.semilogx.LineWidth=1.5;

% Standard semilogy plot
Sd.PlotParam.semilogy.LineWidth=1.5;

% Patch plot
Sd.PlotParam.patch.FaceColor='interp';
Sd.PlotParam.patch.EdgeColor='none';
Sd.PlotParam.patch.LineWidth=1;

% Surf plot
Sd.PlotParam.surf.MeshStyle='both';

% Mesh plot
Sd.PlotParam.mesh.MeshStyle='both';

% Quiver plot
Sd.PlotParam.quiver.LineWidth=1;

% Set default axis parameters
Sd.AxisParam.Colormap=jet(256);
Sd.AxisParam.Title.String='Solution';

% Set default legend parameters
Sd.LegendParam.Title.String={};
Sd.LegendParam.Visible='off';
Sd.LegendParam.Location='southeast';

% Set default colorbar parameters
Sd.ColorbarParam.Visible='off';

Sd.step_count=true;
Sd.step=1;

[F]=merge(Fd, F);

%% Initialize

n=length(varargin);
S=cell(n,1);

for k=1:n
    if isrow(varargin{k})
        varargin{k}=varargin{k}';
    end

    S{k}=arrayfun(@(x) checkInput(struct_fields,Sd,x), varargin{k});
end

S=cat(1,S{:});
n=length(S);

%% Post-processing

nb_steps=arrayfun(@(x) numel(x.DynamicParam), S);
nb_steps_min=min(nb_steps);
nb_steps_max=max(nb_steps);

if F.shrink
    map=round(1+(nb_steps-1)./(nb_steps_min-1).*(0:(nb_steps_min-1)));
    target=nb_steps_min;
else
    map=round(1+(nb_steps-1)./(nb_steps_max-1).*(0:(nb_steps_max-1)));
    target=nb_steps_max;
end

% Add fields to struct array
nb_steps=num2cell(nb_steps);
map=num2cell(map, 2);

[S.nb_steps]=nb_steps{:};
[S.map]=map{:};

% Checking time step
check_step=[S.step]>[S.nb_steps];

if any(check_step)
    warning('Step chosen exceeds number of steps. Changing step to maximum step.')
    max_step=num2cell([S(check_step).nb_steps]);
    [S(check_step).step]=max_step{:};
end

% Add legend to figure
legends=[S.LegendParam];

if ismember('on', {legends.Visible})
    F.LegendParam.Visible='on';
end

%% Plotting interface

Plot=cell(n,1);
Axis=cell(n,1);
Colorbar=cell(n,1);

fig=figure;
fig=setObject(fig, F.FigParam);

t=tiledlayout(1,n);
t=setObject(t, F.TiledParam);

for k=1:n
    axis=nexttile;
    axis.NextPlot='add'; % Hold on

    switch S(k).type
        
        case {'plot', 'loglog', 'semilogx', 'semilogy'}
            Param=concatenate(S(k).PlotParam(1), S(k).DynamicParam(1), {'XData', 'YData'}); % XData and YData must be specified before anything else
            p=feval(S(k).type, axis, Param.XData, Param.YData);
            p=setObject(p, Param, {'XData', 'YData'});

        case 'patch'
            p=patch(axis, S(k).PlotParam(1), S(k).DynamicParam(1));

        case 'mesh'
            p=mesh(axis, S(k).DynamicParam(1).ZData);
            p=setObject(p, S(k).PlotParam(1));

        case 'surf'
            id=structfun(@iscell, S(k).DynamicParam(1));
            l=structfun(@length, S(k).DynamicParam(1));

            if any(id)
                np=unique(l(id)); % Number of plots
                for j=1:np
                    p(j)=feval(S(k).type, axis);
                end
            else
                p=feval(S(k).type, axis);
            end

            p=setObject(p, S(k).PlotParam(1));
            p=setObject(p, S(k).DynamicParam(1));

        case 'quiver'
            p=quiver(axis, S(k).PlotParam(1), S(k).DynamicParam(1));

    end

    c=colorbar(axis);
    l=legend(axis, S(k).LegendParam.Title.String); % Title must be specified before any other parameter
    axis=setObject(axis, S(k).AxisParam);
    c=setObject(c, S(k).ColorbarParam);
    l=setObject(l, S(k).LegendParam, {'Title'});
    refreshdata(axis);

    Plot{k}=p;
    Axis{k}=axis;
    Colorbar{k}=c;
end

% Set common legend for entire figure
lg=legend([Plot{1}], F.LegendParam.Title.String);
lg=setObject(lg, F.LegendParam, {'Title'});

if F.view_video % View the solution over the entire simulation

    frames(target)=struct('cdata',[],'colormap',[]);
    frames(1)=getframe(fig);
    for step = 1:target
        pause(0.01)
        for k=1:n
            Title={S(k).AxisParam.Title.String, [' at step ' num2str(S(k).map(step)) ' out of ' num2str(S(k).map(target))]};
            Plot{k}=setObject(Plot{k}, S(k).DynamicParam(S(k).map(step)));
            Axis{k}.Title.String=[Title{[true, S(k).step_count]}];
            refreshdata(Axis{k});
        end
        frames(step)=getframe(fig);
    end

    if F.write_video
        video=VideoWriter(F.vid_name, 'MPEG-4');
        video.FrameRate=(target)/F.vid_duration;
        video.Quality=100;
        open(video);
        writeVideo(video,frames);
        close(video);
    end

else % View the solution at a particular time snapshot
    for k=1:n
        Title={S(k).AxisParam.Title.String, [' at step ' num2str(S(k).step) ' out of ' num2str(S(k).nb_steps)]};
        Plot{k}=setObject(Plot{k}, S(k).DynamicParam(S(k).step));
        Axis{k}.Title.String=[Title{[true, S(k).nb_steps>1 && S(k).step_count]}];
        refreshdata(Axis{k});
        box(Axis{k}, 'on')
    end
end


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

                    switch class(param.(fields{i}))
                        case 'cell' % Objects take different parameters
                            [object.(fields{i})]=param.(fields{i}){:};
                        case 'double'
                            if size(param.(fields{i}),2)==s % Objects take different parameters (stored in a regular array)
                                a=num2cell(param.(fields{i}),1);
                                [object.(fields{i})]=a{:};
                            else % All objects take the same parameter
                                [object.(fields{i})]=deal(param.(fields{i}));
                            end
                        otherwise % All objects take the same parameter
                            [object.(fields{i})]=deal(param.(fields{i}));
                    end
                end
            end
        end

    end

%% Check Input

    function[So]=checkInput(fields, Sd, Si)

        if any(~ismember(fieldnames(Si), fields))
            error('Supplied fields are not supported')
        end

        So=merge(Sd, Si, {'PlotParam', 'DynamicParam'});

        if isfield(Si, 'PlotParam')
            So.PlotParam=merge(Sd.PlotParam.(So.type), Si.PlotParam);
        else
            So.PlotParam=Sd.PlotParam.(So.type);
        end

        if isfield(Si, 'DynamicParam')
            So.DynamicParam=Si.DynamicParam;
        else
            So.DynamicParam=So.PlotParam;
        end

    end

%% Concatenate and sort array

    function[Ss]=concatenate(S1,S2,names)

        % Merge the struct arrays S1 and S2 and sort the fields following
        % the order in the cell array "names".

        Ss=merge(S1,S2);
        fields=fieldnames(Ss);
        ofields=[names'; fields];
        ofields=unique(ofields, 'stable');
        Ss=orderfields(Ss,ofields);

    end
end