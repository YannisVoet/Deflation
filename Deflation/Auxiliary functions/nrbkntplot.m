function nrbkntplot (nurbs, options)

% NRBKNTPLOT: Plot a NURBS entity with the knots subdivision.
% 
% Calling Sequence:
% 
%   NRBKNTPLOT(nurbs)
%   NRBKNTPLOT(nurbs, name, value)
% 
% INPUT:
% 
%   nurbs:      NURBS curve, surface or volume, see nrbmak.
%   nsub:       Number of evaluation points, for a surface or volume, a row 
%               vector with the number of points along each direction. 
%               Default: 50.
%   colormap:   Colormap (string or [R G B] vector). Default: 'summer'.
%   light:      Lighting for plot ('on' or 'off'). Default: 'on'.
%
% Example:
%
%   Plot the test surface with its knot vector
%
%   NRBKNTPLOT(nrbtestsrf)
%
% See also:
% 
%   nrbctrlplot
%
%    Copyright (C) 2011, 2012, 2021 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

arguments
    nurbs
    options.light {mustBeText} = 'on'
    options.colormap = 'summer'
    options.nsub {mustBeNumeric} = 50
end

if (nargin < 1)
  error ('nrbkntplot: Need a NURBS to plot!');
end

light=options.light;
cmap=options.colormap;
nsub=options.nsub;

colormap (cmap);

hold_flag = ishold;

if (iscell (nurbs.knots))
 if (size (nurbs.knots,2) == 2) % plot a NURBS surface
   if (nargin < 2)
     nsub = [50 50];
   elseif (numel(nsub) == 1)
     nsub = [nsub nsub];
   end
   nrbplot (nurbs, nsub, 'light', light, 'colormap', cmap);
   hold on

   % And plot the knots
   knt1 = unique (nurbs.knots{1}(nurbs.order(1):end-nurbs.order(1)+1));
   knt2 = unique (nurbs.knots{2}(nurbs.order(2):end-nurbs.order(2)+1));
   p1 = nrbeval (nurbs, {knt1, linspace(knt2(1),knt2(end),nsub(2)+1)});
   p2 = nrbeval (nurbs, {linspace(knt1(1),knt1(end),nsub(1)+1), knt2});

  if (any (nurbs.coefs(3,:)))
    % surface in a 3D space
    for ii = 1:numel(knt1)
      plot3 (squeeze(p1(1,ii,:)), squeeze(p1(2,ii,:)), squeeze(p1(3,ii,:)),'k');
    end
    for ii = 1:numel(knt2)
      plot3 (squeeze(p2(1,:,ii)), squeeze(p2(2,:,ii)), squeeze(p2(3,:,ii)),'k'); 
    end
  else
    % plain surface
    for ii = 1:numel(knt1)
      plot (squeeze(p1(1,ii,:)), squeeze (p1(2,ii,:)),'k'); 
    end
    for ii = 1:numel(knt2)
      plot (p2(1,:,ii),p2(2,:,ii),'k');
    end
  end



 elseif (size (nurbs.knots,2) == 3) % plot a NURBS volume
   if (nargin < 2)
     nsub = [25 25 25];
   elseif (numel(nsub) == 1)
     nsub = [nsub nsub nsub];
   end
   % Plot the boundaries
   bnd = nrbextract (nurbs);
   nrbkntplot (bnd(1), nsub(2:3));
   hold on
   for iface = 2:6
     inds = setdiff(1:3, ceil(iface/2));
     nrbkntplot (bnd(iface), nsub(inds));
   end
 end
else % plot a NURBS curve
  if (nargin < 2)
    nsub = 1000;
  end
  nrbplot (nurbs, nsub);
  hold on

  % And plot the knots
   order = nurbs.order;
   p = nrbeval (nurbs, unique (nurbs.knots(order:end-order+1)));

   if (any (nurbs.coefs(3,:))) % plot a 3D curve
     plot3 (p(1,:), p(2,:), p(3,:), 'rx'); 
   else                     % plot a 2D curve
     plot (p(1,:), p(2,:), 'rx'); 
   end

end

if (~hold_flag)
  hold off
end

end

