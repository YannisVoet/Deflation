function[M]=sppagemtimes(A,B)

% SPPAGEMTIMES: This generalizes Matlab's pagemtimes functions when one of
% the matrices is sparse.
%
% Z = SPPAGEMTIMES(X,Y) computes the pagewise multiplication when one of
% the arrays is a sparse matrix.

if issparse(A)
    M=arrayfun(@(i) A*B(:,:,i), 1:size(B,3), 'UniformOutput', false);
else
    M=arrayfun(@(i) A(:,:,i)*B, 1:size(A,3), 'UniformOutput', false);
end

M=cat(3,M{:});
end