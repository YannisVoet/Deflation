% MSH_SPLIT: Splits the approximation space into large and small subspaces.
%
% CALLING SEQUENCE:
%
% [IL, IS] = MSH_SPLIT(sp_trimmed, msh_trimmed) returns the index sets for
% large and small basis functions. The resulting large and small subspaces
% are defined, respectively, as VL = span(B(IL)) and VS = span(B(IS)).
%
% [IL, IS, ILl, ISl] = MSH_SPLIT(sp_trimmed, msh_trimmed) returns the index 
% sets in the local numbering.
%
% [IL, IS, ILl, ISl, good_ids, bad_ids] = MSH_SPLIT(sp_trimmed, msh_trimmed) 
% also returns the indices of large and small elements.
%
% [IL, IS] = MSH_SPLIT(sp_trimmed, msh_trimmed, gamma) provides a parameter
% gamma in (0,1] for partitioning the mesh into large and small elements.
% An element is considered small if |T ∩ Omega| < theta |T|. Default: 0.1
% If gamma = 1, all cut elements are considered small.
%
% References:
% [1] A. Buffa, R. Puppi, and R. Vázquez. A minimal stabilization procedure
% for isogeometric methods on trimmed geometries. SINUM, 2020.
% [2] E. Burman, P. Hansbo, M. G. Larson, and K. Larsson. Extension operators
% for trimmed spline spaces. CMAME, 2023.
%
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

function[ILg, ISg, ILl, ISl, good_ids, bad_trimmed_ids] = msh_split(sp_trimmed, msh_trimmed, gamma)

arguments
    sp_trimmed
    msh_trimmed
    gamma {mustBeGreaterThan(gamma,0), mustBeLessThanOrEqual(gamma,1)} = 0.1
end

%% Algorithm

% Parition the mesh into large and small elements
trimmed_ids = msh_trimmed.reparam.trim_elem_ids;
mesh_trimmed = msh_evaluate_element_list(msh_trimmed.msh_cart, trimmed_ids);
el_size_trimmed = (sum (mesh_trimmed.quad_weights .* abs (mesh_trimmed.jacdet), 1));

vol_trimmed = zeros (msh_trimmed.reparam.nb_trim_elems, 1);
for iel = 1:msh_trimmed.reparam.nb_trim_elems
    for jel = 1 : msh_trimmed.reparam.trim_elems(iel).nb_tiles
        mesh_tile = msh_tile(msh_trimmed, iel, jel);
        vol_trimmed(iel) = vol_trimmed(iel) + sum (mesh_tile.jacdet .* mesh_tile.quad_weights, 1);
    end
end

good_trimmed_ids = trimmed_ids (vol_trimmed ./ el_size_trimmed(:) >= gamma);
bad_trimmed_ids = setdiff(trimmed_ids, good_trimmed_ids);
good_ids = union(good_trimmed_ids, msh_trimmed.reparam.non_trim_elem_ids);

% Indices in global numbering
ILg=sp_get_basis_functions(sp_trimmed.space_untrimmed, msh_trimmed.msh_cart, good_ids);
ISg=setdiff(sp_trimmed.active_dofs, ILg);

% Indices in local numbering
ILl=sp_trimmed.global_to_active(ILg);
ISl=sp_trimmed.global_to_active(ISg);

end