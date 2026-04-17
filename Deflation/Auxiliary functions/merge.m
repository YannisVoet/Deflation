function[Param]=merge(Param,param,varargin)

% MERGE: merges the fields of struct arrays.
%
% S = MERGE(S1,S2) merges the fields of struct arrays S1 and S2. If S1
% and S2 contain fields in common, those in S1 are overriden by those in
% S2.
%
% S = MERGE(S1,S2,X) excludes the fields contained in the cell array X.

if nargin>2
    X=varargin{1};
else
    X={};
end

if ~isempty(param)

    fields=setdiff(fieldnames(param), X, 'stable');

    for k=1:length(fields)

        if isstruct(param.(fields{k}))
            if isfield(Param, fields(k))
                Param.(fields{k})=merge(Param.(fields{k}), param.(fields{k}));
            else
                Param.(fields{k})=param.(fields{k});
            end
        else
            Param.(fields{k})=param.(fields{k});
        end
    end
end