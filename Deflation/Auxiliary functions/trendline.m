function[c]=trendline(x,y,options)

% TRENDLINE: Plots a trendline for the data.
%
% [c]=TRENDLINE(x,y) returns the coefficient vector for the trendline
% c(1)*x+c(2).
%
% [c]=TRENDLINE(x,y,name,value) specifies optional name/value pair arguments.
% Among them are:
%   'model'     - specifies the model to fit. Available models are:
%       'standard'    - the trendline is c(1)*x+c(2).
%       'loglog'      - the trendline is exp(c(2))*(x.^c(1)).
%                 Default: loglog
%   'linestyle' - line style for the plot (see Matlab documentation).
%                 Default: '--'
%   'linewidth' - line width for the plot (see Matlab documentation).
%                 Default: 2
%   'color'     - line color for the plot (see Matlab documentation).
%                 Default: 'k'
%   'variable'  - variable for the legend. The printed legend is:
%                 'O(round(c(1))*variable' for standard plots.
%                 'O(variable^round(c(1)))' for loglog plots.
%                 Example: 'epsilon'
%                 Default: 'h'
%   'digits'    - rounding for decimal digits.
%                 Default: 0.
%   'subset'    - Removes points for improving the trendline.
%                 Default: false.
%   'tol'       - Euclidean norm error tolerance for point removal.
%                 Default: 1e-1.
%   'maxiter'   - Maximum number of points removed.
%                 Default: min([10 length(x)]).


arguments
    x
    y
    options.model = 'loglog'
    options.linestyle = '--';
    options.linewidth = 2;
    options.color = 'k';
    options.variable = 'h'
    options.digits = 0;
    options.subset=false;
    options.tol=1e-1;
    options.maxiter=min([10 length(x)]);
end

switch options.model
    case 'standard'
        mapx=@(x) x; invmapx=@(x) x;
        mapy=@(x) x; invmapy=@(x) x;

    case 'semilogx'
        mapx=@(x) log(x); invmapx=@(x) exp(x);
        mapy=@(x) x;      invmapy=@(x) x;

    case 'semilogy'
        mapx=@(x) x;      invmapx=@(x) x;
        mapy=@(x) log(x); invmapy=@(x) exp(x);

    case 'loglog'
        mapx=@(x) log(x); invmapx=@(x) exp(x);
        mapy=@(x) log(x); invmapy=@(x) exp(x);

end

% Remove infinite or nan values
idx=~isfinite(y) | isnan(y);

x(idx)=[];
y(idx)=[];

if isempty(x)
    return
end

% Turn to column vectors
x=reshape(x, [], 1);
y=reshape(y, [], 1);
n=length(x);

if options.subset

    X=mapx(x);
    Y=mapy(y);

    c=polyfit(X, Y, 1)'; % Least squares solution
    err=vecnorm([X ones(n,1)]*c-Y, 2, 2);
    niter=1;
    actid=1:n;
    I=[];

    while max(err)>options.tol && niter<=options.maxiter

        [~,Imax]=max(err);

        X(Imax)=[]; % Remove element with largest error
        Y(Imax)=[];

        n=length(X);
        c=polyfit(X, Y, 1)'; % Least squares solution
        err=vecnorm([X ones(n,1)]*c-Y, 2, 2);
        niter=niter+1;
        I=[I actid(Imax)];
        actid(Imax)=[];
    end

    x(I)=[];
    y(I)=[];
end

c=polyfit(mapx(x), mapy(y), 1)';

switch options.model
    case 'standard'
        plot(x, 0.9*c(2)+c(1)*x, 'LineStyle', options.linestyle, 'Color', options.color, 'LineWidth', options.linewidth, 'DisplayName', ['$O(' num2str(round(c(1),options.digits)) options.variable ')$'])
    case 'semilogx'
        semilogx(x, 0.99*c(2)+c(1)*mapx(x), 'LineStyle', options.linestyle, 'Color', options.color, 'LineWidth', options.linewidth, 'DisplayName', ['$O(' num2str(round(c(1),options.digits)) options.variable ')$'])
    case 'semilogy'
        semilogy(x, exp(c(2)-1e-1)*(x.^c(1)), 'LineStyle', options.linestyle, 'Color', options.color, 'LineWidth', options.linewidth, 'DisplayName', ['$O(' options.variable '^{' num2str(round(c(1),options.digits)) '})$'])
    case 'loglog'
        loglog(x, exp(c(2)-1e-1)*(x.^c(1)), 'LineStyle', options.linestyle, 'Color', options.color, 'LineWidth', options.linewidth, 'DisplayName', ['$O(' options.variable '^{' num2str(round(c(1),options.digits)) '})$'])
end



