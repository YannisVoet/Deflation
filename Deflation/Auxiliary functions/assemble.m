function[A]=assemble(V, mapi, mapj, ni, nj)

% ASSEMBLE: assembles element (or single-patch) matrices into a single
% global matrix.
%
% A = ASSEMBLE(V,map) assembles the local matrices contained in V into a
% global matrix A. The element matrices may either be stored in vectorized
% format such that V = [vect(E1),...,vect(Em)] or as elements of a cell
% array such that V = {E1,...,Em}. The vectorized format is typically used
% when all element matrices have the same size (e.g. classical FEM) whereas
% the cell array format allows for element matrices of different sizes
% (e.g. multi-patch IGA). The map provided maps local degrees of freedom to
% global ones. Its format must be consistent with the format of V:
% If V is in vectorized format, map is a connectivity matrix with as many
% columns as elements. If V is a cell array, map = {map1,...,mapm} is also
% a cell array whose elements are vectors.
%
% A = ASSEMBLE(V,mapi,mapj) assembles the local matrices for different
% connectivity matrices for rows and columns.
%
% A = ASSEMBLE(V,mapi,mapj,ni,nj) also specifies the size of A as an
% ni x nj matrix.

arguments
    V
    mapi
    mapj = mapi
    ni = 'auto'
    nj = ni
end


switch class(V)
    case 'double'   % V = [vect(E1),...,vect(Em)] (all Ej have the same size)

        mapi=check_size(V, mapi);
        mapj=check_size(V, mapj);

        nne=sqrt(size(V,1));
        e=ones(nne,1);
        I=kron(e, mapi);
        J=kron(mapj,e);

        if isequal(ni, 'auto')
            ni=max(mapi, [], 'all');
        end

        if isequal(nj, 'auto')
            nj=max(mapj, [], 'all');
        end

        A=sparse(I(:), J(:), V(:), ni, nj);

    case 'cell'     % V = {E1,...,Em} (Ej could have different sizes)

        m=length(V);

        rs=cell(m,1);
        cs=cell(m,1);
        vs=cell(m,1);

        for j=1:m
            [rs{j}, cs{j}, vs{j}]=find(V{j});
            rs{j}=mapi{j}(rs{j});
            cs{j}=mapj{j}(cs{j});
        end

        if isequal(ni, 'auto')
            ni=max(cat(1,rs{:}));
        end

        if isequal(nj, 'auto')
            nj=max(cat(1,cs{:}));
        end

        A=sparse(cat(1,rs{:}), cat(1,cs{:}), cat(1,vs{:}), ni, nj);
end
end


function[map]=check_size(V,map)

if size(map,2)~=size(V,2)
    if size(map,1)==size(V,2)
        map=map';
    else
        error('The connectivity matrix must contain as many columns as elements.')
    end
end

end



% [A1,A2,...,An] = ASSEMBLE(V1,map1,V2,map2...,Vn,mapn) assembles Vk into a
% global matrix Ak according to mapk for k = 1,2,...,n.

%% Old version
% n=length(varargin);
%
% if mod(n,2)==0
%     n=n/2;
% else
%     error('Input must be provided in pairs.')
% end
%
% varargout=cell(1,n);
%
% for k=1:n
%
%     V=varargin{2*k-1};
%     map=varargin{2*k};
%
%     switch class(V)
%         case 'double'   % V = [vect(E1),...,vect(Em)] (all Ej have the same size)
%
%             if size(map,2)~=size(V,2)
%                 if size(map,1)==size(V,2)
%                     map=map';
%                 else
%                     error('The connectivity matrix must contain as many columns as elements.')
%                 end
%             end
%
%             nne=sqrt(size(V,1));
%             e=ones(nne,1);
%             I=kron(e, map);
%             J=kron(map,e);
%
%             varargout{k}=sparse(I(:), J(:), V(:));
%
%         case 'cell'     % V = {E1,...,Em} (Ej could have different sizes)
%
%             m=length(V);
%
%             rs=cell(m,1);
%             cs=cell(m,1);
%             vs=cell(m,1);
%
%             for j=1:m
%                 [rs{j}, cs{j}, vs{j}]=find(V{j});
%                 rs{j}=map{j}(rs{j});
%                 cs{j}=map{j}(cs{j});
%             end
%
%             varargout{k}=sparse(cat(1,rs{:}), cat(1,cs{:}), cat(1,vs{:}));
%     end
% end