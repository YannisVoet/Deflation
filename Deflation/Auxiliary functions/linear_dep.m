function ISr = linear_dep(sp_trimmed, msh_trimmed, gamma, threshold)

% LINEAR_DEP: Detect nearly linearly dependent basis functions by support
% overlap.
%
% ISr = LINEAR_DEP(sp_trimmed, msh_trimmed) returns the indices
% (in local/active numbering) of basis functions to deflate, using a
% support overlap criterion.
%
% ISr = LINEAR_DEP(sp_trimmed, msh_trimmed, gamma) only deflates the
% elements whose volume fraction is smaller than gamma in (0,1].
% Default: 1. Note: changing this parameter is not recommended.
%
% ISr = LINEAR_DEP(sp_trimmed, msh_trimmed, gamma, threshold) specifies
% the overlap threshold. Two functions are flagged if the following
% conditions are met:
%  1) The functions overlap.
%  2) The volume fraction of the cut elements in their support is
% (on average) smaller than the threshold. Default threshold: 0.25.

arguments
    sp_trimmed
    msh_trimmed
    gamma {mustBeGreaterThan(gamma,0), mustBeLessThanOrEqual(gamma,1)} = 1
    threshold {mustBeGreaterThan(threshold,0), mustBeLessThanOrEqual(threshold,1)} = 0.25
end

ISr = [];
if (isempty(msh_trimmed.reparam.trim_elems))
    return;
end

trimmed_ids = msh_trimmed.reparam.trim_elem_ids;
nb_trim_elems = msh_trimmed.reparam.nb_trim_elems;
mesh_trimmed = msh_evaluate_element_list(msh_trimmed.msh_cart, trimmed_ids);
el_size_trimmed = sum(mesh_trimmed.quad_weights .* abs(mesh_trimmed.jacdet), 1);

ndof = sp_trimmed.space_untrimmed.ndof;
nsh_max = sp_trimmed.space_untrimmed.nsh_max;
V = [];

vol_trimmed = zeros(nb_trim_elems, 1);
for iel = 1:nb_trim_elems
    for jel = 1:msh_trimmed.reparam.trim_elems(iel).nb_tiles
        mesh_tile = msh_tile(msh_trimmed, iel, jel);
        vol_trimmed(iel) = vol_trimmed(iel) + sum(mesh_tile.jacdet .* mesh_tile.quad_weights, 1);
    end
    func_ids = sp_get_basis_functions(sp_trimmed.space_untrimmed, msh_trimmed.msh_cart, trimmed_ids(iel));
    V = [V func_ids];
end

% Numerator always uses active (trimmed) volumes.
active_weight = vol_trimmed(:)';

% Connectivity matrix S: rows = global basis functions, columns = trimmed elements.
elems = repmat(1:nb_trim_elems, nsh_max, 1);
S = sparse(V(:), elems(:), ones(nsh_max * nb_trim_elems, 1), ndof, nb_trim_elems);

% Keep only small functions returned by msh_split.
[~, IS, ~, ISl] = msh_split(sp_trimmed, msh_trimmed, gamma);
S = S(IS, :);

% Remove rows with no support on trimmed elements.
nonzero_rows = any(S, 2);
S = S(nonzero_rows, :);
ISl = ISl(nonzero_rows);


if isempty(ISl)
    return;
end

S_binary = S > 0;


% Within each group, check pairs for support overlap above threshold.
flagged = [];
n=size(S_binary,1);

for ii = 1:n
    supp_i = S_binary(ii, :);
    for jj = (ii+1):n
        supp_j = S_binary(jj, :);
        if any((supp_i & supp_j))
            shared = (supp_i | supp_j) * active_weight';
            denom = (supp_i | supp_j) * el_size_trimmed';
            if shared > 0 && shared / denom <= threshold
                flagged = [flagged [ISl(ii); ISl(jj)]];
            end
        end
    end
end


if ~isempty(flagged)
    ISr = unique(flagged(:));
end

end
