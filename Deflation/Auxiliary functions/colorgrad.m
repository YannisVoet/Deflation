function[c]=colorgrad(c1,c2,n)

% COLORGRAD: Returns a color gradient for any intermdiate color between a
% light color c1 and a dark color c2.
%
% c = COLORGRAD(c1,c2) returns the function handle c(t) = t*c2+(1-t)*c1
% where t is in [0, 1] and c1,c2 are color triplets.
% Hence, c(0) = c1 and c(1) = c2.
%
% c = COLORGRAD(c1,c2,n) returns instead a sampling of c(t) at n points ti 
% uniformly distributed between 0 and 1 such that c = [c(t1);...;c(tn)]

c = @(t) t*c2+(1-t)*c1;

if nargin > 2
    s=linspace(0,1,n)';
    c=c(s);
end
end