function O = tricomm_apply_inits(O, items)
% TRICOMM_APPLY_INITS  Fold "--init NAME=VALUE" arguments into an options struct.
%
%   O = tricomm_apply_inits(O, {'Ar=20', 'Bs=6', 'Ab=0.36'})
%
% NAME is a compartment of the 10-state model:
%
%   As Ar   sensitive / resistant LM
%   Bs Br   sensitive / resistant PM
%   Cs Cr   sensitive / resistant PL
%   P L M   free amino acids
%   Ab      antibiotic (an alias for Ab_init)
%
% Values are in 1e7 CFU/mL for the strains. Each one overrides whatever the
% community preset and the experimental data would have supplied for that
% compartment, and nothing else -- so seeding a resistant sub-population does
% not disturb the rest of the initial state.

NAMES = {'As','Ar','Bs','Br','Cs','Cr','P','L','M'};

for i = 1:numel(items)
    it = strtrim(items{i});
    eq = find(it == '=', 1);
    if isempty(eq)
        error('--init expects NAME=VALUE, got "%s".', it);
    end
    name = strtrim(it(1:eq-1));
    x    = str2double(strtrim(it(eq+1:end)));
    if ~isfinite(x)
        error('--init %s: "%s" is not a finite number.', name, strtrim(it(eq+1:end)));
    end
    if x < 0
        error('--init %s: a negative initial level (%g) is not meaningful.', name, x);
    end

    if strcmp(name, 'Ab')
        O.Ab_init = x;
    elseif any(strcmp(name, NAMES))
        O.(['init_' name]) = x;
    else
        error('--init: unknown compartment "%s".\n  Valid: %s Ab', ...
              name, strjoin(NAMES, ' '));
    end
end
end
