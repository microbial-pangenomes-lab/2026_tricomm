function P = tricomm_apply_sets(P, sets)
% TRICOMM_APPLY_SETS  Apply "--set key=value" overrides to a parameter struct.
%
%   P = tricomm_apply_sets(P, {'alpha_C=0.06', 'delta=*1.1'})
%
% A value starting with '*' is a MULTIPLIER of the current value, so
% "--set alpha_A=*0.5" halves it. Overriding a key the parameter file never
% defines is an error, for the same reason unknown keys are rejected there.

for i = 1:numel(sets)
    s  = sets{i};
    eq = find(s == '=', 1);
    if isempty(eq)
        error('--set expects key=value, got "%s".', s);
    end
    key = strtrim(s(1:eq-1));
    val = strtrim(s(eq+1:end));

    if ~isfield(P, key)
        error('--set: unknown parameter "%s".', key);
    end

    if ~isempty(val) && val(1) == '*'
        m = str2double(val(2:end));
        if ~isfinite(m)
            error('--set %s: "%s" is not a valid multiplier.', key, val);
        end
        P.(key) = P.(key) * m;
    else
        x = str2double(val);
        if ~isfinite(x)
            error('--set %s: "%s" is not a finite number.', key, val);
        end
        P.(key) = x;
    end
end
end
