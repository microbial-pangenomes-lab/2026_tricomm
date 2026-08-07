function v = tricomm_opt(opt, key, dflt, mode)
% TRICOMM_OPT  Read one option out of a tricomm_argparse struct.
%
%   v = tricomm_opt(opt, 'tmax', 1100, 'num')   last --tmax as a number
%   v = tricomm_opt(opt, 'out', 'x.csv')        last --out as a string
%   v = tricomm_opt(opt, 'set', {}, 'all')      every --set, as a cellstr
%   v = tricomm_opt(opt, 'quiet', false, 'bool')
%
% mode: 'str' (default), 'num', 'bool', 'all'. Absent options return dflt.

if nargin < 4, mode = 'str'; end

if ~isfield(opt, key)
    v = dflt;
    return
end
vals = opt.(key);

switch mode
    case 'all'
        v = vals;
    case 'bool'
        s = lower(vals{end});
        v = ~any(strcmp(s, {'false', '0', 'no', 'off'}));
    case 'num'
        v = str2double(vals{end});
        if ~isfinite(v)
            error('option --%s: "%s" is not a number.', key, vals{end});
        end
    otherwise
        v = vals{end};
end
end
