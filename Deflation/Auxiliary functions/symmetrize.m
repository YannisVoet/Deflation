function[varargout] = symmetrize(varargin)

% SYMMETRIZE: Returns the symmetric part of a matrix M.
%
% Ms = SYMMETRIZE(M) returns the symmetric part Ms of M.
%
% [Ms1,...,Msn] = SYMMETRIZE(M1,...,Mn) returns the symmetric part Msi of Mi 
% for i = 1,...,n.

n=nargin;
varargout=cell(1,n);

for k=1:n
    varargout{k} = .5*(varargin{k} + varargin{k}');
    assert(norm(varargin{k} - varargout{k}, 'fro') / norm(varargin{k}, 'fro') < 1e-10)
end
end

