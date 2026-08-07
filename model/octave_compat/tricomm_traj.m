function M = tricomm_traj(S, vars)
% TRICOMM_TRAJ  Pull named state variables out of a tricomm_simulate result.
%
%   M     = tricomm_traj(S, {'LM','PM','PL'})
%   names = tricomm_traj()          % the valid names, for validating input
%
% Returns an nt x numel(vars) matrix. Recognised names:
%
%   LM PM PL      strain totals, sensitive + resistant   (aliases A, B, C)
%   As Ar         sensitive / resistant LM
%   Bs Br         sensitive / resistant PM
%   Cs Cr         sensitive / resistant PL
%   P L M         free amino acids
%   Ab            antibiotic
%   total         LM + PM + PL
%
% All strain values are in units of 1e7 CFU/mL, as everywhere else in the model.
% The raw-state names are the columns of S.Y in the order dYdt_model.m returns
% them; the totals are what the published figures actually plot.

RAW = {'As','Ar','Bs','Br','Cs','Cr','P','L','M','Ab'};

if nargin < 2
    % Name list only, so callers can reject a typo before spending a sweep on it.
    M = [{'LM','PM','PL','total','A','B','C'}, RAW];
    return
end

if ischar(vars), vars = {vars}; end
M = zeros(numel(S.t), numel(vars));

for i = 1:numel(vars)
    v = strtrim(vars{i});
    switch v
        case {'LM','A'},  M(:,i) = S.A;
        case {'PM','B'},  M(:,i) = S.B;
        case {'PL','C'},  M(:,i) = S.C;
        case 'total',     M(:,i) = S.A + S.B + S.C;
        otherwise
            k = find(strcmp(v, RAW), 1);
            if isempty(k)
                error(['tricomm_traj: unknown variable "%s".\n' ...
                       '  Valid: LM PM PL total %s (and aliases A B C).'], ...
                      v, strjoin(RAW, ' '));
            end
            M(:,i) = S.Y(:,k);
    end
end
end
