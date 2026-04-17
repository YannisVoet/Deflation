function s = findspan_new(U,u,n)                
% FINDSPAN  Find the span of a B-Spline knot vector at a parametric point
%
% Calling Sequence:
% 
%   s = findspan(U,u,n)
% 
%  INPUT:
% 
%    U - knot sequence
%    u - parametric points
%    n - number of control points (or basis functions)
%    
% 
%  OUTPUT:
% 
%    s - knot span index
%
%  Modification of Algorithm A2.1 from 'The NURBS BOOK' pg68
%
%    Copyright (C) 2010 Rafael Vazquez
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

if max(u(:))>U(end) || min(u(:))<U(1)
  error('Some value is outside the knot span')
end

% Making sure U is a row vector
if size(U,1)>1
    U=U';
end

% Making sure u is a colum vector
if size(u,2)>1
    u=u';
end


[~,s]=find(cumsum(U<=u,2,'reverse')==1);
s(u==U(n+1))=n; % For open knot vectors, the maximum value of mu is n (the number of basis functions)
end