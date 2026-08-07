function [pos, opt] = tricomm_argparse(args, boolflags)
% TRICOMM_ARGPARSE  Minimal "--flag value" parser for the run_* CLI scripts.
%
%   [pos, opt] = tricomm_argparse(argv(), {'help','quiet'})
%
% Accepts  --key value  and  --key=value. Names in boolflags take no value.
% Dashes inside a flag become underscores, so --solve-alpha-c reads back as
% opt.solve_alpha_c.
%
% pos is a cellstr of positional arguments. opt is a struct whose every field is
% a CELLSTR holding each occurrence in order, so repeatable flags (--set,
% --sweep) work without special-casing. Read it with tricomm_opt.

if nargin < 2, boolflags = {}; end
boolflags = strrep(boolflags, '-', '_');

pos = {};
opt = struct();
i   = 1;

while i <= numel(args)
    a = args{i};

    if numel(a) < 2 || ~strcmp(a(1:2), '--')
        pos{end+1} = a;
        i = i + 1;
        continue
    end

    body = a(3:end);
    eq   = find(body == '=', 1);
    if ~isempty(eq)
        key = body(1:eq-1);
        val = body(eq+1:end);
        i   = i + 1;
    else
        key = body;
        if any(strcmp(strrep(key, '-', '_'), boolflags))
            val = 'true';
            i   = i + 1;
        else
            if i + 1 > numel(args)
                error('option --%s needs a value.', key);
            end
            val = args{i+1};
            i   = i + 2;
        end
    end

    key = strrep(key, '-', '_');
    if isempty(key)
        error('malformed option "%s".', a);
    end
    if isfield(opt, key)
        opt.(key){end+1} = val;
    else
        opt.(key) = {val};
    end
end
end
