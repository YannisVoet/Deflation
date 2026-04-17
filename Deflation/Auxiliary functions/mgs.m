function[Q, varargout]=mgs(V, varargin)

% MGS: Modified Gram-Schmidt algorithm. Based on Algorithm 5.2.6 of [1].
%
% Q = MGS(V) computes an orthonornal basis for the columns of the 
% (full rank) matrix V. 
%
% Q = MGS(V,A) computes an A-orthonormal basis instead. A is a positive
% definite matrix either provided explicitly or as a function handle.
%
% [Q,R] = MGS(V,...) also returns an upper triangular matrix such that 
% V = Q*R.
%
% Reference:
% [1] G. H. Golub and C. F. Van Loan. Matrix computations. JHU press, 2013.

%% Set algorithm parameters

[n,m]=size(V);

Param = inputParser;
Param.addRequired('V');
Param.addOptional('A', @(x) x, @(x) isa(x, 'double') || isa(x, 'function_handle'));
Param.parse(V, varargin{:});

%% Retrieve parameters
A=Param.Results.A;

if isa(A, 'double')
    A=@(x) A*x;
end

%% Modified Gram-Schmidt

R=zeros(m,m);
Q=zeros(n,m);

for k=1:m
    R(k,k)=sqrt(V(:,k)'*A(V(:,k)));
    Q(:,k)=V(:,k)/R(k,k);

    for j=k+1:m
        R(k,j)=Q(:,k)'*A(V(:,j));
        V(:,j)=V(:,j)-Q(:,k)*R(k,j);
    end
end
varargout{1}=R;
end