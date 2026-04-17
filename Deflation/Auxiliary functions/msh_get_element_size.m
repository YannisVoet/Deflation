% MSH_GET_ELEMENT_SIZE: Returns the size of trimmed elements.
%
% USAGE:
%
%  [vol_trim, size_trim, ids_trim, size_non_trim, ids_non_trim] = msh_get_element_size(msh)
%
% INPUT:
%
% msh:   mesh object (see msh_trimming)
%       
% OUTPUT:
%
% vol_trim:         size (or volume) of the trimmed elements intersected with the
%                   physical domain; i.e. |T ∩ Omega|.
% size_trim:        size (or volume) of the trimmed elements; i.e. |T|
% ids_trim:         trimmed element ids.
% size_non_trim:    size (or volume) of the non-trimmed elements.
% ids_non_trim:     non-trimmed element ids
%
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


function [vol_trimmed, el_size_trim, trim_elem_ids, el_size_non_trim, non_trim_elem_ids] = msh_get_element_size(msh_trimmed)

% Non-trimmed elements
non_trim_elem_ids = msh_trimmed.reparam.non_trim_elem_ids;
non_trim_elems = msh_evaluate_element_list(msh_trimmed.msh_cart, non_trim_elem_ids);
el_size_non_trim = (sum (non_trim_elems.quad_weights .* abs (non_trim_elems.jacdet), 1))';

% Trimmed elements
trim_elem_ids = msh_trimmed.reparam.trim_elem_ids;
trim_elems = msh_evaluate_element_list(msh_trimmed.msh_cart, trim_elem_ids);
el_size_trim = (sum (trim_elems.quad_weights .* abs (trim_elems.jacdet), 1))';

vol_trimmed = zeros (msh_trimmed.reparam.nb_trim_elems, 1);
for iel = 1:msh_trimmed.reparam.nb_trim_elems
    for jel = 1 : msh_trimmed.reparam.trim_elems(iel).nb_tiles
        mesh_tile = msh_tile(msh_trimmed, iel, jel);
        vol_trimmed(iel) = vol_trimmed(iel) + sum (mesh_tile.jacdet .* mesh_tile.quad_weights, 1);
    end
end

end