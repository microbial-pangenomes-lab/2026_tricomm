%% run_model.m  --  Run the tricomm model with a chosen parameter set and
%%                   save the full time course to a CSV/TSV file.
%
% This is a plotting-free, Octave-friendly version of
% 1_simulation/linear_scale_version/tricomm_lin.m. Edit the CONFIG block
% below, then run it (see octave_compat/README.md for the exact command).
%
% Output: one row per time point, columns:
%   time, As, Ar, Bs, Br, Cs, Cr, P, L, M, Ab, A_tot, B_tot, C_tot,
%   A_CFUmL, B_CFUmL, C_CFUmL
% where *_tot = sensitive + resistant (model units) and *_CFUmL = *_tot*1e7.
%
% Two ways to set parameters:
%   1. Edit the defaults in the CONFIG block below and run `run_model`.
%   2. Pre-set any of them, then run, without editing this file, e.g.:
%        octave-cli -p octave_compat --eval \
%          "community=0; TMAX=1000; OUTFILE='nodrug.csv'; run('octave_compat/run_model.m')"
%      Every CONFIG value uses the pre-set variable if it exists, else the default.

% ======================= CONFIG (defaults) =======================
% Each line: keep a pre-set value if the caller defined one, else use default.
if ~exist('DATA_DIR','var'),     DATA_DIR   = '../data';           end  % path to data/ folder
if ~exist('OUTFILE','var'),      OUTFILE    = 'model_output.csv';  end  % .csv->comma, .tsv->tab

% --- Experiment / community selection (see tricomm_lin.m header for meaning) ---
if ~exist('antibiotic','var'),   antibiotic = 1;   end  % 1 = ampicillin 16xMIC, 2 = colistin
if ~exist('dose','var'),         dose       = 2;   end  % only for antibiotic==2 (1=1xMIC,2=2xMIC)
if ~exist('community','var'),    community  = 3;   end  % 0 = no-drug, 1..4 = drug communities

% --- Model parameters (the "set of parameters" you want to run) ---
if ~exist('epsilon','var'),      epsilon      = 1;      end  % antibiotic-induced killing rate
if ~exist('eta','var'),          eta          = 0.005;  end  % conjugation (plasmid transfer) rate
if ~exist('Ab_init','var'),      Ab_init      = 0.18;   end  % initial antibiotic concentration
if ~exist('scale_factor','var'), scale_factor = 1.1;    end  % resistant growth multiplier (plasmid cost)

if ~exist('mu_A','var'),    mu_A    = 0.0003;  end
if ~exist('mu_B','var'),    mu_B    = 0.0005;  end
if ~exist('mu_C','var'),    mu_C    = 0.0004;  end
if ~exist('alpha_A','var'), alpha_A = 0.7645;  end
if ~exist('alpha_B','var'), alpha_B = 0.0567;  end
if ~exist('alpha_C','var'), alpha_C = 0.1202;  end
if ~exist('delta','var'),   delta   = 0.0692;  end

if ~exist('TMAX','var'),    TMAX    = [];       end  % final time (h). []->last experimental point.
% =================================================================

%% ---- Load experimental data (for initial conditions) ----
cfu_AMP = readmatrix(fullfile(DATA_DIR, 'pOXA48_interpolated.csv'));
cfu_COL = readmatrix(fullfile(DATA_DIR, 'PN23_interpolated.csv'));

switch antibiotic
    case 1, time = cfu_AMP(:,1);
    case 2, time = cfu_COL(:,1);
end

switch community
    case 0
        cfu_A = cfu_AMP(:,2);  cfu_B = cfu_AMP(:,3);  cfu_C = cfu_AMP(:,4);
    case 1
        if antibiotic == 1
            cfu_A = cfu_AMP(:,5);  cfu_B = cfu_AMP(:,6);  cfu_C = cfu_AMP(:,7);
        else
            if dose == 1, cols = [2 3 4]; else, cols = [5 6 7]; end
            cfu_A = cfu_COL(:,cols(1)); cfu_B = cfu_COL(:,cols(2)); cfu_C = cfu_COL(:,cols(3));
        end
    case 2
        if antibiotic == 1
            cfu_A = cfu_AMP(:,8);  cfu_B = cfu_AMP(:,9);  cfu_C = cfu_AMP(:,10);
        else
            cfu_A = cfu_COL(:,8);  cfu_B = cfu_COL(:,9);  cfu_C = cfu_COL(:,10);
        end
    case 3
        if antibiotic == 1
            cfu_A = cfu_AMP(:,11); cfu_B = cfu_AMP(:,12); cfu_C = cfu_AMP(:,13);
        else
            if dose == 1, cols = [11 12 13]; else, cols = [14 15 16]; end
            cfu_A = cfu_COL(:,cols(1)); cfu_B = cfu_COL(:,cols(2)); cfu_C = cfu_COL(:,cols(3));
        end
    case 4
        cfu_A = cfu_AMP(:,14); cfu_B = cfu_AMP(:,15); cfu_C = cfu_AMP(:,16);
end

%% ---- Assemble parameters (with community-specific overrides) ----
mu_Ar = mu_A * scale_factor;
mu_Br = mu_B * scale_factor;
mu_Cr = mu_C * scale_factor;

if community == 0
    epsilon = 0; eta = 0; Ab_init = 0;
end
if community == 1 || community == 4
    eta = 0;
end

params = [mu_A; mu_Ar; mu_B; mu_Br; mu_C; mu_Cr; ...
          alpha_A; alpha_B; alpha_C; delta; epsilon; eta];

%% ---- Initial conditions ----
switch community
    case {0, 1}
        As_init = cfu_A(1);  Ar_init = 0;
    case {2, 3, 4}
        As_init = 0;         Ar_init = cfu_A(1);
end
Y0 = [As_init; Ar_init; cfu_B(1); 0; cfu_C(1); 0; 0; 0; 0; Ab_init];

if isempty(TMAX), TMAX = time(end); end

%% ---- Solve ----
[t, Y] = ode45(@(t,Y) dYdt_model(t, Y, params, community, antibiotic), 1:TMAX, Y0);

As = Y(:,1); Ar = Y(:,2); Bs = Y(:,3); Br = Y(:,4); Cs = Y(:,5); Cr = Y(:,6);
P  = Y(:,7); L  = Y(:,8); M  = Y(:,9); Ab = Y(:,10);
A_tot = As + Ar;  B_tot = Bs + Br;  C_tot = Cs + Cr;

%% ---- Write output table ----
header = {'time','As','Ar','Bs','Br','Cs','Cr','P','L','M','Ab', ...
          'A_tot','B_tot','C_tot','A_CFUmL','B_CFUmL','C_CFUmL'};
out = [t, As, Ar, Bs, Br, Cs, Cr, P, L, M, Ab, ...
       A_tot, B_tot, C_tot, A_tot*1e7, B_tot*1e7, C_tot*1e7];

[~, ~, ext] = fileparts(OUTFILE);
if strcmpi(ext, '.tsv'), delim = '\t'; else, delim = ','; end

fid = fopen(OUTFILE, 'w');
if fid < 0, error('Cannot open %s for writing.', OUTFILE); end
fprintf(fid, [strjoin(header, delim) '\n']);
fmt = [repmat(['%.10g' delim], 1, size(out,2)-1) '%.10g\n'];
fprintf(fid, fmt, out.');   % transpose: fprintf walks column-major
fclose(fid);

fprintf('Wrote %d rows x %d cols to %s\n', size(out,1), size(out,2), OUTFILE);
fprintf('  antibiotic=%d community=%d  epsilon=%g eta=%g Ab_init=%g scale=%g\n', ...
        antibiotic, community, epsilon, eta, Ab_init, scale_factor);
fprintf('  final CFU/mL:  A(LM)=%.4g  B(PM)=%.4g  C(PL)=%.4g\n', ...
        A_tot(end)*1e7, B_tot(end)*1e7, C_tot(end)*1e7);
