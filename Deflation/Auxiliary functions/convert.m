function[varargout]=convert(A, format, varargin)

% CONVERT: Converts storage format of factor matrices.
%
% A = CONVERT(A,fout) converts the array A to the format specified by
% fout. Supported input and output formats are 'tensor', 'cell', 'vect'.
% Format description:
% 'tensor'  - the factor matrices are stored along the pages of a third
%             order tensor.
% 'cell'    - the factor matrices are stored as elements of cell arrays.
% 'vect'    - the factor matrices are stored in vectorized format as
%             [vect(A1),...,vect(Ar)]. If A is provided in vectorized
%             format, it is assumed that all Ak are n x n matrices for
%             k = 1,...,r and r < n^2.
%
% [A,s] = CONVERT(A,fout) also returns the size of A as s = [n,m,r], where
% [n,m] is the size of the factor matrices and r is the number of factor
% matrices.
%
% [A,s,fin] = CONVERT(A,fout) also returns the input format of the array A.
%
% [A1,A2,...,Am,...] = CONVERT(A,{f1,f2,...,fm}) converts the array A to an
% array Aj with format fj for j = 1,...,m.

if ~iscell(format)
    format={format};
end

% First vectorize all
switch class(A)
    case 'double'
        s=size(A);

        if ismatrix(A)
            if s(1)==s(2)
                f='tensor';
                s=[s 1];
            else
                f='vect';
                s=[sqrt(s(1)) sqrt(s(1)) s(2)];
            end
        else
            f='tensor';
        end

        A=reshape(A, [prod(s(1:2)) s(3)]);

    case 'cell'
        r=length(A);
        sizes=cellfun(@size, A, 'UniformOutput', false);
        sizes=cat(1,sizes{:});

        n=unique(sizes(:,1));
        m=unique(sizes(:,2));

        if  isscalar(n) && isscalar(m)
            s=[n m r];
            f='cell';
            A=cell2vect(A);
        else
            error('Conversion process is impossible. Check sizes.')
        end
end

% Then convert vectorized arrays to desired format
nf=length(format);
varargout=cell(nf+2, 1);

for k=1:nf
    switch format{k}
        case 'vect'
            varargout{k}=A;
        case 'tensor'
            if s(3)==1
                varargout{k}=reshape(A, s);
            else
                varargout{k}=reshape(full(A), s);
            end
        case 'cell'
            varargout{k}=mat2cell(reshape(A, [s(1) prod(s(2:3))]), s(1), s(2)*ones(1,s(3)));
    end
end

varargout{nf+1}=s;
varargout{nf+2}=f;