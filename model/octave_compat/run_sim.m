%% run_sim.m -- run the TriComm model for one parameter file and save the run.
%
%   octave-cli model/octave_compat/run_sim.m PARFILE [options]
%
% See model/README.md, or run with --help.

here = fileparts(mfilename('fullpath'));
addpath(here);

BOOLS = {'help', 'quiet', 'mse', 'solve_alpha_c', 'report_only', ...
         'drug', 'no_drug', 'conjugation', 'no_conjugation'};
[pos, opt] = tricomm_argparse(argv(), BOOLS);

HELP = {
'run_sim.m -- integrate the TriComm model for one parameter file.'
''
'  octave-cli model/octave_compat/run_sim.m PARFILE [options]'
''
'Options'
'  --out FILE         output path; .csv -> comma, .tsv -> tab'
'                     (default <label>_c<community>.csv in the current folder)'
'  --community N      preset: 0 = no drug (default), 1..4 = drug communities.'
'                     Sets the drug and conjugation switches, whether LM starts'
'                     sensitive or resistant, and which data columns seed the'
'                     run. Every part of it can be overridden below.'
'  --antibiotic N     1 = ampicillin 16xMIC, 2 = colistin'
'  --dose N           1 = 1xMIC, 2 = 2xMIC (antibiotic 2 only)'
'  --tmax H           final time in hours (default: last experimental point)'
'  --tstart H         first output time (default 1)'
'  --dt H             output step (default 1)'
'  --solver NAME      ode45 (default, matches the published scripts) or lsode'
'  --extinct-below X  a strain whose TOTAL (sensitive+resistant) falls below X'
'                     CFU/mL is set to zero for the rest of the run, so'
'                     extinction is absorbing. Default 0 = off (published'
'                     behaviour). Use 1 for one cell per mL.'
'  --rtol X / --atol X  solver tolerances. Unset = the solver default, which'
'                     for ode45 is what the published runs used. ode45 defaults'
'                     to AbsTol 1e-6 = 10 CFU/mL here, so tighten both when the'
'                     question turns on small populations.'
''
'Antibiotic'
'  --drug / --no-drug   force it on or off, whatever the community says'
'  --ab0 X            initial antibiotic level (default 0.18)'
'  --ab-in X          antibiotic in the feed. Defaults to --ab0, i.e. the drug'
'                     is replenished as in the experiments and settles at'
'                     delta*ab0/(gamma+delta), just below the dose. Pass'
'                     --ab-in 0 for the published single pulse.'
'  --epsilon X        killing rate per unit antibiotic (default 1)'
'  --gamma X          antibiotic decay rate (default 0.0032 ampicillin,'
'                     0.0050 colistin)'
''
'Plasmid'
'  --conjugation / --no-conjugation   force transfer on or off'
'  --eta X            conjugation rate (default 0.005)'
''
'Initial levels'
'  --y0 A,B,C         strain totals in 1e7 CFU/mL, overriding the data. The'
'                     community decides whether A lands in As or Ar.'
'  --init NAME=VALUE  one compartment at a time, overriding both the preset and'
'                     --y0 (repeatable). NAME is As Ar Bs Br Cs Cr P L M or Ab,'
'                     e.g. --init Ar=2 --init Bs=6 seeds a resistant LM'
'                     sub-population without touching anything else.'
''
'  --set key=value    override any parameter; value may be *m for a multiplier'
'                     (repeatable, e.g. --set alpha_C=0.06 --set delta=*1.1).'
'                     epsilon, eta, Ab_init, gamma and init_* are parameters'
'                     too, so run_sweep.m can sweep them.'
'  --solve-alpha-c    replace alpha_C by the root of a*b*c = a+b+c+2, i.e. put'
'                     the set exactly on the neutral surface R = 1'
'  --mse              also report the fit to the KAN no-drug replicates'
'  --report-only      print the diagnostics, write no file'
'  --quiet            print nothing but errors'
'  --help             this text'
''
'Examples'
'  octave-cli model/octave_compat/run_sim.m model/params/solved.par \'
'      --tmax 1100 --out nodrug_1100h.csv'
'  octave-cli model/octave_compat/run_sim.m model/params/calc_par.par \'
'      --solve-alpha-c --mse --report-only'
''
'  # the no-drug community, but dose it anyway'
'  octave-cli model/octave_compat/run_sim.m model/params/ori_par.par \'
'      --drug --ab0 0.36 --tmax 200'
''
'  # community 3 with the plasmid unable to spread'
'  octave-cli model/octave_compat/run_sim.m model/params/ori_par.par \'
'      --community 3 --no-conjugation --tmax 200'
''
'  # start from a hand-built community: 1%% of LM resistant, no drug'
'  octave-cli model/octave_compat/run_sim.m model/params/ori_par.par \'
'      --no-drug --init As=11.6 --init Ar=0.15 --init Bs=8.7 --init Cs=8'
};

if tricomm_opt(opt, 'help', false, 'bool') || isempty(pos)
    printf('%s\n', HELP{:});
    if isempty(pos) && ~tricomm_opt(opt, 'help', false, 'bool')
        error('no parameter file given.');
    end
    return
end
if numel(pos) > 1
    error('expected one parameter file, got %d: %s', numel(pos), strjoin(pos, ' '));
end

quiet = tricomm_opt(opt, 'quiet', false, 'bool');

%% ---- Parameters ---------------------------------------------------------
P = tricomm_load_params(pos{1});
P = tricomm_apply_sets(P, tricomm_opt(opt, 'set', {}, 'all'));

if tricomm_opt(opt, 'solve_alpha_c', false, 'bool')
    [~, aC] = tricomm_R(P);
    if ~isfinite(aC)
        error(['--solve-alpha-c: no positive root exists for this set ' ...
               '(needs alpha_A*alpha_B > delta^2).']);
    end
    P.alpha_C = aC;
end

%% ---- Run ----------------------------------------------------------------
O = struct();
for f = {'community', 'antibiotic', 'dose', 'tmax', 'tstart', 'dt', ...
         'epsilon', 'eta', 'gamma', 'rtol', 'atol', 'extinct_below'}
    if isfield(opt, f{1}), O.(f{1}) = tricomm_opt(opt, f{1}, NaN, 'num'); end
end
if isfield(opt, 'ab0'),   O.Ab_init = tricomm_opt(opt, 'ab0',   NaN, 'num'); end
if isfield(opt, 'ab_in'), O.Ab_in   = tricomm_opt(opt, 'ab_in', NaN, 'num'); end
if isfield(opt, 'solver'), O.solver = tricomm_opt(opt, 'solver', 'ode45'); end
if isfield(opt, 'y0')
    O.y0 = sscanf(strrep(tricomm_opt(opt, 'y0', ''), ',', ' '), '%f').';
end

% --drug / --no-drug and --conjugation / --no-conjugation. Giving both halves of
% a pair is a contradiction rather than a last-one-wins, so say so.
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

S = tricomm_simulate(P, O);
[R, aC_star, frac] = tricomm_R(P);

%% ---- Report -------------------------------------------------------------
if ~quiet
    printf('parameter set : %s   (%s)\n', P.label, P.source);
    printf('  mu    A %-12.6g B %-12.6g C %-12.6g\n', P.mu_A, P.mu_B, P.mu_C);
    printf('  alpha A %-12.6g B %-12.6g C %-12.6g   delta %.6g\n', ...
           P.alpha_A, P.alpha_B, P.alpha_C, P.delta);
    printf('  R = a*b*c/(a+b+c+2) = %.8g   ->  %s\n', R, tricomm_fate(R));
    if isfinite(aC_star)
        printf('  alpha_C for R = 1 exactly: %.6g  (this set is %+.2f%% off)\n', ...
               aC_star, 100 * (P.alpha_C - aC_star) / aC_star);
    end
    printf('  equilibrium composition LM/PM/PL: %.3f / %.3f / %.3f\n', frac);
    printf('run           : community %d, antibiotic %d, dose %d, solver %s\n', ...
           S.community, S.antibiotic, S.dose, S.solver);
    if S.drug
        printf('  drug          ON   Ab_init %.4g, epsilon %.4g, gamma %.4g\n', ...
               S.Ab_init, S.epsilon, S.gamma);
        if S.Ab_in > 0
            printf('  replenished   Ab_in %.4g -> plateau %.4g (%.1f%% of Ab_init)\n', ...
                   S.Ab_in, S.Ab_star, 100 * S.Ab_star / S.Ab_init);
        else
            printf('  replenished   no -- single pulse, half-life %.2f h (--ab-in 0)\n', ...
                   log(2) / (S.gamma + P.delta));
        end
    else
        printf('  drug          OFF\n');
    end
    if S.conjugation
        printf('  conjugation   ON   eta %.4g\n', S.eta);
    else
        printf('  conjugation   OFF\n');
    end
    printf('  t = %g .. %g h, %d points\n', S.t(1), S.t(end), numel(S.t));
    printf('  initial state As %.4g Ar %.4g | Bs %.4g Br %.4g | Cs %.4g Cr %.4g\n', ...
           S.y0(1), S.y0(2), S.y0(3), S.y0(4), S.y0(5), S.y0(6));
    printf('  start CFU/mL  LM %.4g  PM %.4g  PL %.4g\n', ...
           (S.y0(1) + S.y0(2)) * 1e7, (S.y0(3) + S.y0(4)) * 1e7, ...
           (S.y0(5) + S.y0(6)) * 1e7);
    printf('  final CFU/mL  LM %.4g  PM %.4g  PL %.4g\n', ...
           S.A(end) * 1e7, S.B(end) * 1e7, S.C(end) * 1e7);
    tot  = S.A + S.B + S.C;
    printf('  total biomass %.6g -> %.6g  (%+.2f%% over the run)\n', ...
           tot(1), tot(end), 100 * (tot(end) - tot(1)) / tot(1));
    % Second half only: by then the transient from the initial condition has
    % died out, so this is the long-run direction rather than the start-up.
    half = max(2, round(numel(S.t) / 2));
    dt_h = S.t(end) - S.t(half);
    if dt_h > 0 && tot(half) > 0 && tot(end) > 0
        printf('  drift over the second half: %+.4g decades / 1000 h\n', ...
               1000 * (log10(tot(end)) - log10(tot(half))) / dt_h);
    end
    if min([S.A; S.B; S.C]) < 1e-6
        printf('  WARNING: at least one strain drops below 1e-6 (1e1 CFU/mL).\n');
    end
end

if tricomm_opt(opt, 'mse', false, 'bool')
    MSE = tricomm_mse(P);
    if ~quiet
        printf('  MSE vs KAN no-drug replicates (<=310 h, log10): %.6g\n', MSE);
    end
end

%% ---- Write --------------------------------------------------------------
if ~tricomm_opt(opt, 'report_only', false, 'bool')
    % Default name comes from the FILE, not the label: labels are prose and
    % contain spaces and brackets.
    [~, base] = fileparts(P.source);
    outfile   = tricomm_opt(opt, 'out', sprintf('%s_c%d.csv', base, S.community));

    header = {'time','As','Ar','Bs','Br','Cs','Cr','P','L','M','Ab', ...
              'A_tot','B_tot','C_tot','A_CFUmL','B_CFUmL','C_CFUmL'};
    out = [S.t, S.Y, S.A, S.B, S.C, S.A * 1e7, S.B * 1e7, S.C * 1e7];

    [~, ~, ext] = fileparts(outfile);
    if strcmpi(ext, '.tsv'), delim = sprintf('\t'); else, delim = ','; end

    fid = fopen(outfile, 'w');
    if fid < 0, error('cannot open %s for writing.', outfile); end
    fprintf(fid, '%s\n', strjoin(header, delim));
    fmt = [repmat(['%.10g' delim], 1, size(out, 2) - 1) '%.10g\n'];
    fprintf(fid, fmt, out.');      % transpose: fprintf walks column-major
    fclose(fid);

    if ~quiet
        printf('wrote %d rows x %d cols to %s\n', size(out, 1), size(out, 2), outfile);
    end
end
