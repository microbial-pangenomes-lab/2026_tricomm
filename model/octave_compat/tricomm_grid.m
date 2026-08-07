function [key, vals] = tricomm_grid(spec, P)
% TRICOMM_GRID  Expand one "--sweep key=spec" argument into a list of values.
%
%   [key, vals] = tricomm_grid('alpha_C=0.05:0.005:0.07', P)
%
% Supported specs, all relative to parameter struct P where marked:
%
%   start:step:stop     inclusive linear range   0.05:0.005:0.07
%   log:start:stop:n    n points, log-spaced     log:1e-4:10:6
%   v1,v2,v3            explicit list            0.05,0.0592,0.0693
%   *m1,*m2,*m3         multipliers of P.key     *0.5,*1,*1.5
%
% Multipliers and literals can be mixed in one list. A single-value spec is
% legal (it just pins the parameter).

eq = find(spec == '=', 1);
if isempty(eq)
    error('--sweep expects key=spec, got "%s".', spec);
end
key  = strtrim(spec(1:eq-1));
body = strtrim(spec(eq+1:end));

if ~isfield(P, key)
    error('--sweep: unknown parameter "%s".', key);
end

if numel(body) > 4 && strcmpi(body(1:4), 'log:')
    parts = strsplit(body(5:end), ':');
    if numel(parts) ~= 3
        error(['--sweep %s: a log range must be log:start:stop:n (three fields ' ...
               'after "log:"), got "%s".'], key, body);
    end
    a = str2double(parts{1});
    b = str2double(parts{2});
    n = str2double(parts{3});
    if ~all(isfinite([a b n]))
        error('--sweep %s: bad log range "%s".', key, body);
    end
    if a <= 0 || b <= 0
        error(['--sweep %s: a log range needs positive endpoints, got %g and ' ...
               '%g.'], key, a, b);
    end
    if n < 1 || n ~= fix(n)
        error('--sweep %s: n must be a positive whole number, got %g.', key, n);
    end
    if n == 1
        vals = a;
    else
        vals = 10 .^ linspace(log10(a), log10(b), n);
        % linspace hits both ends, but the powers of 10 come back as 0.1000...1
        % and similar; snap them so the output columns read cleanly.
        vals([1 end]) = [a, b];
    end
elseif any(body == ':')
    parts = strsplit(body, ':');
    if numel(parts) ~= 3
        error(['--sweep %s: a range must be start:step:stop (three fields), ' ...
               'got "%s".'], key, body);
    end
    a = str2double(parts{1});
    s = str2double(parts{2});
    b = str2double(parts{3});
    if ~all(isfinite([a s b])) || s == 0
        error('--sweep %s: bad range "%s".', key, body);
    end
    if sign(b - a) ~= sign(s) && b ~= a
        error('--sweep %s: step %g never reaches %g from %g.', key, s, b, a);
    end
    vals = a:s:b;
    % Include the endpoint when the step does not divide the interval exactly.
    if ~isempty(vals) && abs(vals(end) - b) > 1e-12 * max(1, abs(b))
        vals(end+1) = b;
    end
else
    items = strsplit(body, ',');
    vals  = zeros(1, numel(items));
    for i = 1:numel(items)
        it = strtrim(items{i});
        if ~isempty(it) && it(1) == '*'
            m = str2double(it(2:end));
            if ~isfinite(m)
                error('--sweep %s: bad multiplier "%s".', key, it);
            end
            vals(i) = P.(key) * m;
        else
            x = str2double(it);
            if ~isfinite(x)
                error('--sweep %s: "%s" is not a finite number.', key, it);
            end
            vals(i) = x;
        end
    end
end

if isempty(vals)
    error('--sweep %s: spec "%s" produced no values.', key, body);
end
end
