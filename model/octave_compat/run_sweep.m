%% run_sweep.m -- sweep parameters over a grid and tabulate what changes.
%
%   octave-cli model/octave_compat/run_sweep.m PARFILE --sweep key=spec [options]
%
% See model/README.md, or run with --help.

here = fileparts(mfilename('fullpath'));
addpath(here);

BOOLS = {'help', 'quiet', 'no_mse', 'no_summary', ...
         'drug', 'no_drug', 'conjugation', 'no_conjugation'};
[pos, opt] = tricomm_argparse(argv(), BOOLS);

HELP = {
'run_sweep.m -- sweep parameters over a grid and tabulate the result.'
''
'  octave-cli model/octave_compat/run_sweep.m PARFILE --sweep key=spec [options]'
''
'Sweep specs (--sweep is repeatable; N of them give an N-D grid)'
'  --sweep alpha_C=0.05:0.005:0.07    start:step:stop, inclusive'
'  --sweep alpha_C=0.0592,0.0693      explicit list'
'  --sweep alpha_A=*0.5,*1,*1.5       multipliers of the value in PARFILE'
'  Any key of the parameter file can be swept: mu_A mu_B mu_C alpha_A alpha_B'
'  alpha_C delta epsilon eta Ab_init Ab_in gamma scale_factor, and the initial'
'  levels'
'  init_As init_Ar init_Bs init_Br init_Cs init_Cr init_P init_L init_M.'
''
'Run length'
'  --tmax H           last time point, h (default 1100)'
'  --tstart H         first time point, h (default 1)'
'  --dt H             output step, h (default 1; raise it to shrink --traj)'
''
'Output'
'  --out FILE         one summary row per run (default sweep_<label>.tsv)'
'  --traj FILE        the full time courses in tidy form: one row per (run,'
'                     time), the swept parameters repeated as columns. This is'
'                     the 3_grid_map view -- every curve it would have plotted.'
'                       run  epsilon  eta  time  Ar  Bs  Br  Cs  Cr'
'  --traj-vars LIST   which variables go in --traj (default LM,PM,PL).'
'                     Also: total As Ar Bs Br Cs Cr P L M Ab'
'  --no-summary       skip the summary file (only makes sense with --traj)'
'  .tsv gives tab-separated, any other extension comma-separated.'
''
'Model options'
'  --community N      preset: 0 = no drug (default), 1..4 = drug communities.'
'                     Sets the drug and conjugation switches and the initial'
'                     state; the flags below override any part of it.'
'  --drug / --no-drug                 force the antibiotic on or off'
'  --conjugation / --no-conjugation   force plasmid transfer on or off'
'  --ab0 X / --epsilon X / --gamma X  antibiotic level, killing rate, decay'
'  --ab-in X          antibiotic in the feed. Defaults to --ab0 (replenished,'
'                     settling at delta*ab0/(gamma+delta)). --ab-in 0 gives the'
'                     published single pulse. Sweeping Ab_init moves both.'
'  --eta X            conjugation rate'
'  --init NAME=VALUE  one initial compartment: As Ar Bs Br Cs Cr P L M Ab'
'                     (repeatable)'
'  --antibiotic N / --dose N          as in run_sim.m'
'  --y0 A,B,C         strain totals in 1e7 CFU/mL, overriding the data'
'  --solver NAME      lsode (default here: ~5x faster) or ode45'
'  --rtol X / --atol X  solver tolerances; unset = the solver default. lsode'
'                     uses 1e-8/1e-11, ode45 its own 1e-3/1e-6 (= 10 CFU/mL).'
'  --set key=value    override a parameter before sweeping (repeatable)'
'  --no-mse           skip the fit to the KAN data (roughly halves the runtime)'
'  --quiet            write the files, print nothing'
'  --help             this text'
''
'Summary columns'
'  run        run id, matching the --traj column suffix'
'  the swept parameter(s), then'
'  R          coexistence ratio a*b*c/(a+b+c+2); >1 grows, <1 washes out'
'  fate       grow / washout / neutral, from R'
'  MSE        mean squared log10 residual vs the KAN no-drug replicates'
'  A_end,B_end,C_end   CFU/mL at tmax'
'  tot_end    total biomass at tmax, 1e7 units'
'  slope      log10 change in total biomass per 1000 h over the second half'
'  minABC     smallest strain abundance reached at any time (extinction check)'
''
'Examples'
'  # the 3_grid_map sweep: 10 x 10 epsilon-eta grid of the drug community,'
'  # every trajectory in one tidy TSV'
'  octave-cli model/octave_compat/run_sweep.m model/params/ori_par.par \'
'      --community 3 --tmax 200 \'
'      --sweep epsilon=0.1:0.1:1 --sweep eta=0.00001:0.00011:0.001 \'
'      --traj grid_map.tsv --traj-vars Ar,Bs,Br,Cs,Cr --no-mse'
''
'  # how sharp is the R = 1 knife edge?'
'  octave-cli model/octave_compat/run_sweep.m model/params/calc_par.par \'
'      --sweep alpha_C=0.055:0.001:0.07'
};

if tricomm_opt(opt, 'help', false, 'bool') || isempty(pos)
    printf('%s\n', HELP{:});
    if isempty(pos) && ~tricomm_opt(opt, 'help', false, 'bool')
        error('no parameter file given.');
    end
    return
end

quiet   = tricomm_opt(opt, 'quiet',  false, 'bool');
do_mse  = ~tricomm_opt(opt, 'no_mse', false, 'bool');
tmax    = tricomm_opt(opt, 'tmax',   1100, 'num');
tstart  = tricomm_opt(opt, 'tstart', 1,    'num');
dt      = tricomm_opt(opt, 'dt',     1,    'num');

P0 = tricomm_load_params(pos{1});
P0 = tricomm_apply_sets(P0, tricomm_opt(opt, 'set', {}, 'all'));

specs = tricomm_opt(opt, 'sweep', {}, 'all');
if isempty(specs)
    error('nothing to sweep: pass at least one --sweep key=spec (see --help).');
end

keys = cell(1, numel(specs));
grid = cell(1, numel(specs));
for i = 1:numel(specs)
    [keys{i}, grid{i}] = tricomm_grid(specs{i}, P0);
end
if numel(unique(keys)) ~= numel(keys)
    error('the same parameter is swept twice: %s.', strjoin(keys, ', '));
end

%% ---- Output selection ---------------------------------------------------
trajfile = tricomm_opt(opt, 'traj', '');
do_traj  = ~isempty(trajfile);
do_sum   = ~tricomm_opt(opt, 'no_summary', false, 'bool');
if ~do_sum && ~do_traj
    error('--no-summary with no --traj would compute everything and write nothing.');
end

tvars = {'LM', 'PM', 'PL'};
if isfield(opt, 'traj_vars')
    tvars = strtrim(strsplit(tricomm_opt(opt, 'traj_vars', ''), ','));
    tvars = tvars(~cellfun(@isempty, tvars));
end
% Check the names now rather than inside the loop, where the failed-run handler
% would turn a typo into a column of NaN instead of an error.
bad = tvars(~ismember(tvars, tricomm_traj()));
if ~isempty(bad)
    error('--traj-vars: unknown variable(s) %s.\n  Valid: %s', ...
          strjoin(bad, ', '), strjoin(tricomm_traj(), ' '));
end

%% ---- Fixed run options --------------------------------------------------
O = struct();
O.tmax   = tmax;
O.tstart = tstart;
O.dt     = dt;
O.solver = tricomm_opt(opt, 'solver', 'lsode');
for f = {'community', 'antibiotic', 'dose', 'epsilon', 'eta', 'gamma', 'rtol', 'atol'}
    if isfield(opt, f{1}), O.(f{1}) = tricomm_opt(opt, f{1}, NaN, 'num'); end
end
if isfield(opt, 'ab0'),   O.Ab_init = tricomm_opt(opt, 'ab0',   NaN, 'num'); end
if isfield(opt, 'ab_in'), O.Ab_in   = tricomm_opt(opt, 'ab_in', NaN, 'num'); end
if isfield(opt, 'y0')
    O.y0 = sscanf(strrep(tricomm_opt(opt, 'y0', ''), ',', ' '), '%f').';
end

if isfield(opt, 'drug') && isfield(opt, 'no_drug')
    error('--drug and --no-drug are contradictory.');
end
if isfield(opt, 'conjugation') && isfield(opt, 'no_conjugation')
    error('--conjugation and --no-conjugation are contradictory.');
end
if isfield(opt, 'drug'),           O.drug = 1; end
if isfield(opt, 'no_drug'),        O.drug = 0; end
if isfield(opt, 'conjugation'),    O.conjugation = 1; end
if isfield(opt, 'no_conjugation'), O.conjugation = 0; end

O = tricomm_apply_inits(O, tricomm_opt(opt, 'init', {}, 'all'));

% A fixed --epsilon / --eta / --init in O would silently outrank the swept value
% from the parameter struct, so refuse the combination rather than sweep a
% parameter that cannot move.
for j = 1:numel(keys)
    ok = keys{j};
    if isfield(O, ok)
        error(['--sweep %s conflicts with the fixed --%s given on the same ' ...
               'command line: the fixed value would win for every run.'], ok, ok);
    end
end

% The same trap one level down: with the drug off, epsilon / Ab_init / gamma are
% forced to zero before they reach the solver, so sweeping them produces N
% identical runs. Same for eta with conjugation off. Catch it here instead of
% letting someone read meaning into a flat sweep.
sweep_community = 0;
if isfield(O, 'community'), sweep_community = O.community;
elseif isfield(P0, 'community'), sweep_community = P0.community; end
drug_on = double(sweep_community ~= 0);
if isfield(O, 'drug'), drug_on = O.drug; elseif isfinite(P0.drug), drug_on = P0.drug; end
conj_on = double(any(sweep_community == [2 3]));
if isfield(O, 'conjugation'), conj_on = O.conjugation;
elseif isfinite(P0.conjugation), conj_on = P0.conjugation; end

for j = 1:numel(keys)
    if ~drug_on && any(strcmp(keys{j}, {'epsilon', 'Ab_init', 'Ab_in', 'gamma'}))
        error(['--sweep %s has no effect: the drug is off (community %g). ' ...
               'Add --drug to turn it on.'], keys{j}, sweep_community);
    end
    if ~conj_on && strcmp(keys{j}, 'eta')
        error(['--sweep eta has no effect: conjugation is off (community %g). ' ...
               'Add --conjugation to turn it on.'], sweep_community);
    end
end

% The output grid tricomm_simulate will use. Every run shares it, which is what
% makes the wide trajectory table well defined.
tvec = (tstart:dt:tmax).';
if tvec(end) ~= tmax, tvec(end+1) = tmax; end
nt = numel(tvec);

%% ---- Build the combination table ----------------------------------------
% First key varies slowest, so the rows read like the nested loops in
% 3_grid_map/grid_map.m.
lens = cellfun(@numel, grid);
n    = prod(lens);
combos = zeros(n, numel(grid));
rep = n;
for j = 1:numel(grid)
    v   = grid{j}(:);
    rep = rep / lens(j);
    combos(:, j) = repmat(kron(v, ones(rep, 1)), n / (rep * lens(j)), 1);
end

if ~quiet
    printf('sweeping %s over %d combination(s)\n', strjoin(keys, ' x '), n);
    printf('  t = %g:%g:%g h (%d points), solver %s, community %g\n', ...
           tstart, dt, tmax, nt, O.solver, tricomm_opt(opt, 'community', 0, 'num'));
    if do_traj
        printf('  trajectories: %s -> %s (%d rows)\n', strjoin(tvars, ','), ...
               trajfile, n * nt);
    end
end

%% ---- Sweep --------------------------------------------------------------
res   = NaN(n, 8);       % R MSE A_end B_end C_end tot_end slope minABC
fates = cell(n, 1);
if do_traj
    TR = NaN(nt, n * numel(tvars));
end

for i = 1:n
    P = P0;
    for j = 1:numel(keys)
        P.(keys{j}) = combos(i, j);
    end

    R = tricomm_R(P);
    [~, fates{i}] = tricomm_fate(R);
    res(i, 1) = R;

    if do_mse
        try
            res(i, 2) = tricomm_mse(P);
        catch
            res(i, 2) = NaN;
        end
    end

    try
        S   = tricomm_simulate(P, O);
        tot = S.A + S.B + S.C;
        half = max(2, round(numel(S.t) / 2));

        res(i, 3) = S.A(end) * 1e7;
        res(i, 4) = S.B(end) * 1e7;
        res(i, 5) = S.C(end) * 1e7;
        res(i, 6) = tot(end);
        % log10 change per 1000 h over the second half: the transient has died
        % out by then, so this measures the long-run direction, not the start.
        dt_h = S.t(end) - S.t(half);
        if dt_h > 0 && tot(half) > 0 && tot(end) > 0
            res(i, 7) = 1000 * (log10(tot(end)) - log10(tot(half))) / dt_h;
        end
        res(i, 8) = min([S.A; S.B; S.C]);

        if do_traj
            if numel(S.t) ~= nt
                error('solver returned %d points, expected %d.', numel(S.t), nt);
            end
            cols = (i-1) * numel(tvars) + (1:numel(tvars));
            TR(:, cols) = tricomm_traj(S, tvars);
        end
    catch err
        % A failed run leaves its summary row and trajectory columns as NaN
        % rather than aborting the whole grid.
        if ~quiet
            printf('  [%d/%d] solver failed (%s)\n', i, n, err.message);
        end
    end

    if ~quiet && (mod(i, 10) == 0 || i == n)
        printf('  %d/%d\n', i, n);
    end
end

%% ---- Write --------------------------------------------------------------
[~, base] = fileparts(P0.source);          % the file name, not the prose label
runid = arrayfun(@(k) sprintf('run%03d', k), (1:n).', 'UniformOutput', false);

if do_sum
    outfile = tricomm_opt(opt, 'out', sprintf('sweep_%s.tsv', base));
    [~, ~, ext] = fileparts(outfile);
    if strcmpi(ext, '.tsv'), d = sprintf('\t'); else, d = ','; end

    header = [{'run'}, keys, ...
              {'R','fate','MSE','A_end','B_end','C_end','tot_end','slope','minABC'}];
    fid = fopen(outfile, 'w');
    if fid < 0, error('cannot open %s for writing.', outfile); end
    fprintf(fid, '%s\n', strjoin(header, d));
    for i = 1:n
        fprintf(fid, '%s%s', runid{i}, d);
        for j = 1:numel(keys)
            fprintf(fid, '%.10g%s', combos(i, j), d);
        end
        fprintf(fid, '%.10g%s%s%s', res(i,1), d, fates{i}, d);
        fprintf(fid, '%.10g%s%.10g%s%.10g%s%.10g%s%.10g%s%.10g%s%.10g\n', ...
                res(i,2), d, res(i,3), d, res(i,4), d, res(i,5), d, ...
                res(i,6), d, res(i,7), d, res(i,8));
    end
    fclose(fid);
end

if do_traj
    [~, ~, ext] = fileparts(trajfile);
    if strcmpi(ext, '.tsv') || isempty(ext), d = sprintf('\t'); else, d = ','; end

    % Tidy layout: one row per (time, run), the swept parameters repeated on
    % every row as ordinary columns. Straight into a groupby or a ggplot facet,
    % and it does not change shape when another --sweep is added.
    fid = fopen(trajfile, 'w');
    if fid < 0, error('cannot open %s for writing.', trajfile); end
    fprintf(fid, '%s\n', strjoin([{'run'}, keys, {'time'}, tvars], d));

    fmt = [repmat(['%.10g' d], 1, numel(tvars) - 1) '%.10g\n'];
    for i = 1:n
        % The parameter block is identical down a run, so build it once.
        pre = runid{i};
        for j = 1:numel(keys)
            pre = sprintf('%s%s%.10g', pre, d, combos(i, j));
        end
        cols = (i-1) * numel(tvars) + (1:numel(tvars));
        for k = 1:nt
            fprintf(fid, '%s%s%.10g%s', pre, d, tvec(k), d);
            fprintf(fid, fmt, TR(k, cols));
        end
    end
    fclose(fid);
end

%% ---- Print --------------------------------------------------------------
if ~quiet
    printf('\n');
    MAXROWS = 40;
    w = '';
    for j = 1:numel(keys), w = [w sprintf('%12s', keys{j})]; end
    printf('%8s%s %10s %8s %9s %11s %11s %10s\n', 'run', w, 'R', 'fate', 'MSE', ...
           'tot_end', 'slope/kh', 'minABC');
    for i = 1:min(n, MAXROWS)
        v = '';
        for j = 1:numel(keys), v = [v sprintf('%12.6g', combos(i,j))]; end
        printf('%8s%s %10.6f %8s %9.4g %11.5g %11.4g %10.3g\n', runid{i}, v, ...
               res(i,1), fates{i}, res(i,2), res(i,6), res(i,7), res(i,8));
    end
    if n > MAXROWS
        printf('  ... %d more rows (see the file)\n', n - MAXROWS);
    end
    printf('\n');
    if do_sum,  printf('wrote %d rows to %s\n', n, outfile); end
    if do_traj, printf('wrote %d trajectory rows (%d runs x %d times) to %s\n', ...
                       n * nt, n, nt, trajfile); end
end
