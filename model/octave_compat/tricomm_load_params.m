function P = tricomm_load_params(parfile)
% TRICOMM_LOAD_PARAMS  Read a TriComm parameter file into a struct.
%
%   P = tricomm_load_params('../params/solved.par')
%
% FILE FORMAT -- one "key = value" per line. Blank lines are ignored; anything
% after a '#' or '%' is a comment. A trailing ';' on a value is tolerated so
% that lines pasted from MATLAB still work.
%
%   mu_A    = 0.00029
%   alpha_A = 0.4054
%   ...
%
% As a shortcut the whole 7-vector can be given on one line in the same order
% the collaborators use in their emails and scripts:
%
%   parCal = [mu_A, mu_B, mu_C, alpha_A, alpha_B, alpha_C, delta];
%
% Later assignments win, so a parCal line followed by "alpha_C = 0.059201"
% overrides just that one entry.
%
% KEYS
%   required : mu_A mu_B mu_C alpha_A alpha_B alpha_C delta
%   drug     : epsilon eta Ab_init Ab_in scale_factor gamma
%              (Ab_in is the feed concentration; unset = follow Ab_init,
%              i.e. the drug is replenished. Set it to 0 for a single pulse)
%   switches : drug conjugation                      (1 = on, 0 = off; unset
%              means "whatever the community preset says")
%   initial  : init_As init_Ar init_Bs init_Br init_Cs init_Cr init_P init_L
%              init_M                                (1e7 CFU/mL; unset means
%              "take it from the experimental data")
%   run      : community antibiotic dose tmax        (defaults for the run;
%              command-line flags of run_sim.m / run_sweep.m override them)
%   label    : free text, echoed in reports
%
% community is a PRESET: it sets drug, conjugation, which compartments are
% seeded, and which data columns seed them. Any of the keys above override it,
% so the preset is a starting point rather than a straitjacket.
%
% Unknown keys are an error rather than a silent no-op: a typo'd parameter name
% would otherwise leave the default in place and quietly change the result.

if nargin < 1 || isempty(parfile)
    error('tricomm_load_params: no parameter file given.');
end
if exist(parfile, 'file') ~= 2
    error('tricomm_load_params: parameter file "%s" not found.', parfile);
end

fid = fopen(parfile, 'r');
if fid < 0
    error('tricomm_load_params: cannot open "%s".', parfile);
end
raw = fread(fid, Inf, '*char').';
fclose(fid);

NUMKEYS = {'mu_A','mu_B','mu_C','alpha_A','alpha_B','alpha_C','delta', ...
           'epsilon','eta','Ab_init','Ab_in','scale_factor','gamma','rtol','atol', ...
           'extinct_below', ...
           'drug','conjugation', ...
           'init_As','init_Ar','init_Bs','init_Br','init_Cs','init_Cr', ...
           'init_P','init_L','init_M', ...
           'community','antibiotic','dose','tmax'};
REQUIRED = {'mu_A','mu_B','mu_C','alpha_A','alpha_B','alpha_C','delta'};
PARCAL_ORDER = {'mu_A','mu_B','mu_C','alpha_A','alpha_B','alpha_C','delta'};

P = struct();
P.label = '';
P.source = parfile;

raw   = strrep(raw, sprintf('\r'), sprintf('\n'));   % tolerate CR / CRLF files
lines = strsplit(raw, sprintf('\n'));

for i = 1:numel(lines)
    ln = regexprep(lines{i}, '[#%].*$', '');
    ln = strtrim(ln);
    if isempty(ln), continue; end

    eq = find(ln == '=', 1);
    if isempty(eq)
        error('%s line %d: expected "key = value", got "%s".', parfile, i, ln);
    end
    key = strtrim(ln(1:eq-1));
    val = strtrim(regexprep(ln(eq+1:end), ';\s*$', ''));

    if strcmpi(key, 'parCal')
        v = sscanf(strrep(strrep(strrep(val, '[', ' '), ']', ' '), ',', ' '), '%f');
        if numel(v) ~= 7
            error('%s line %d: parCal needs 7 numbers (%s), got %d.', ...
                  parfile, i, strjoin(PARCAL_ORDER, ', '), numel(v));
        end
        for j = 1:7
            P.(PARCAL_ORDER{j}) = v(j);
        end
        continue
    end

    if strcmpi(key, 'label')
        P.label = val;
        continue
    end

    if ~any(strcmp(key, NUMKEYS))
        error('%s line %d: unknown key "%s".\n  Valid keys: %s, parCal, label.', ...
              parfile, i, key, strjoin(NUMKEYS, ', '));
    end

    x = str2double(val);
    if ~isfinite(x)
        error('%s line %d: "%s" is not a finite number for key "%s".', ...
              parfile, i, val, key);
    end
    P.(key) = x;
end

missing = REQUIRED(~cellfun(@(k) isfield(P, k), REQUIRED));
if ~isempty(missing)
    error('%s: missing required parameter(s): %s.', parfile, strjoin(missing, ', '));
end

% Defaults for anything the file did not set. The drug-community values match
% 1_simulation/linear_scale_version/tricomm_lin.m.
% NaN means "not set here": tricomm_simulate then falls back to the community
% preset (drug, conjugation, gamma) or to the experimental data (init_*).
DEFAULTS = {'epsilon', 1; 'eta', 0.005; 'Ab_init', 0.18; 'Ab_in', NaN; ...
            'scale_factor', 1.1; ...
            'gamma', NaN; 'drug', NaN; 'conjugation', NaN; ...
            'rtol', NaN; 'atol', NaN; 'extinct_below', 0; ...
            'init_As', NaN; 'init_Ar', NaN; 'init_Bs', NaN; 'init_Br', NaN; ...
            'init_Cs', NaN; 'init_Cr', NaN; ...
            'init_P', NaN; 'init_L', NaN; 'init_M', NaN; ...
            'community', 0; 'antibiotic', 1; 'dose', 2; 'tmax', NaN};
for i = 1:size(DEFAULTS, 1)
    if ~isfield(P, DEFAULTS{i,1})
        P.(DEFAULTS{i,1}) = DEFAULTS{i,2};
    end
end

if isempty(P.label)
    [~, base] = fileparts(parfile);
    P.label = base;
end
end
