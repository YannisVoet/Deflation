function[M]=kron2mat(varargin)

% KRON2MAT: Explicitly computes the Kronecker product matrix from the
% factor matrices.
%
% M = KRON2MAT(M1,M2,...,Md) computes the Kronecker product matrix from the
% factor matrices stored in Mk; i.e. M = ∑ M1j ⨂ M2j ⨂ ... ⨂ Mdj.
% The arrays Mk may be either third order tensors or cell arrays but must
% contain the same number of factor matrices.

d=length(varargin);
blocksize=cell(1,d);
np=zeros(1,d);

for k=1:d
    if isa(varargin{k}, 'double')
        varargin{k}=convert(varargin{k}, 'cell');
    end

    [m,n]=cellfun(@size, varargin{k}, 'UniformOutput', true);

    if any(m-m(1)) || any(n-n(1))
        error('Inconsistent size of factor matrices.')
    else
        blocksize{k}=[m(1) n(1)];
    end

    np(k)=length(varargin{k});
end

if any(np-np(1))
    error('Inconsistent number of factor matrices.')
else
    r=np(1);
end

varargin=cat(1,varargin{:});
F=arrayfun(@(k) krong(varargin{:,k}), 1:r, 'UniformOutput', false);
M=sumg(F{:});

    function[M]=krong(varargin)

        % KRONG: Kronecker product of multiple factor matrices
        % (generalizing Matlab's kron function).
        %
        % M = KRONG(M1,M2,...,Md) returns M = M1 ⨂ M2 ⨂ ... ⨂ Md

        if nargin==1
            M=varargin{1};
        else
            M=kron(varargin{1}, krong(varargin{2:end}));
        end

    end

    function[M]=sumg(varargin)

        % SUMG: Sums together multiple factor matrices.
        %
        % M = SUMG(M1,M2,...,Md) returns M = M1 + M2 + ... + Md

        if nargin==1
            M=varargin{1};
        else
            M=varargin{1} + sumg(varargin{2:end});
        end

    end

end