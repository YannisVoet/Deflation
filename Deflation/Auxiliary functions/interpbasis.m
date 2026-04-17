function[B,varargout]=interpbasis(space,options)

% INTERPBASIS: Constructs a (truncated) interpolatory spline basis.
%
% B = INTERPBASIS(space) returns the collocation matrix for the spline
% space at the Demko points. See options for other interpolation points.
%
% [B, C] = INTERPBASIS(space) also returns the invese of B. For truncated
% interpolatory splines, all entries of C that are smaller in magnitude to
% a tolerance are substituted with zeros.
%
% [B, C, pts] = INTERPBASIS(space) also returns the interpolation points.
%
% B = INTERPBASIS(space, name, value) specifies name/value pair arguments.
% Available options:
%   'interp_pts'    - Interpolation points, either provided explicitly as a
%                     cell array or chosen among the following list:
%                     'greville', 'demko', 'uniform'. Default: 'demko'.
%   'trunc'         - Basis truncation (true/false). Default: true.
%   'tol'           - Tolerance for basis truncation. Default: 1e-14.

arguments
    space
    options.interp_pts = 'demko'
    options.trunc = true;
    options.tol = 1e-14;
end

% Knot vector
knots=space.knots;
degree=space.degree;
ndof_dir=space.ndof_dir;

r=length(knots);

B=cell(1,r);
C=cell(1,r);
pts=cell(1,r);

for k=1:r

    if iscell(options.interp_pts)
        pts{k}=options.interp_pts{k};
    else
        switch options.interp_pts
            case 'greville'
                pts{k}=aveknt(knots{k}, degree(k)+1);
            case 'demko'
                pts{k}=chbpnt(knots{k}, degree(k)+1);
            case 'uniform'
                pts{k}=linspace(0,1,ndof_dir(k));
        end
    end

    B{k}=spcol(knots{k}, degree(k)+1, pts{k});
    C{k}=inv(B{k});

    if options.trunc
        C{k}(abs(C{k})<options.tol)=0; % Univariate spline truncation
    end

end

B=fliplr(B);
C=fliplr(C);

B=kron2mat(B{:});
C=kron2mat(C{:});

% Convert to sparse matrices
[IB,JB,VB]=find(B);
[IC,JC,VC]=find(C);

B=sparse(IB,JB,VB,space.ndof,space.ndof);
C=sparse(IC,JC,VC,space.ndof,space.ndof);

varargout{1}=C;
varargout{2}=pts;
end