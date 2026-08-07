function S = tricomm_simulate(P, O)
% TRICOMM_SIMULATE  Integrate the 10-state TriComm model for one parameter set.
%
%   S = tricomm_simulate(P)
%   S = tricomm_simulate(P, O)
%
% P is a struct from tricomm_load_params. O is an options struct; every field is
% optional and falls back to P and then to the defaults below:
%
%   community   0 = no drug, 1..4 = drug communities   (default from P)
%   antibiotic  1 = ampicillin 16xMIC, 2 = colistin    (default from P)
%   dose        1 = 1xMIC, 2 = 2xMIC (antibiotic 2)    (default from P)
%   tstart      first output time, h                   (default 1)
%   tmax        last output time, h                    (default: last data point)
%   dt          output step, h                         (default 1)
%   solver      'ode45' or 'lsode'                     (default 'ode45')
%   y0          [A B C] strain totals in 1e7 CFU/mL, overriding the data
%   data_dir    folder holding *_interpolated.csv      (default ../../data/model)
%
% COMMUNITY IS A PRESET, NOT A CONSTRAINT. It decides four things, each of
% which can be overridden independently through P or O:
%
%   drug         antibiotic present?   on unless community 0
%   conjugation  plasmid transfers?    on only for communities 2 and 3
%   the split    is LM seeded sensitive (As) or resistant (Ar)?
%   the data     which columns of the *_interpolated.csv seed the run
%
% Overrides, all optional:
%
%   drug         1 / 0                 force the antibiotic on or off
%   conjugation  1 / 0                 force plasmid transfer on or off
%   epsilon      killing rate per unit Ab      used when drug is on
%   Ab_init      initial antibiotic level      used when drug is on
%   Ab_in        antibiotic in the feed        defaults to Ab_init, so a dose is
%                one number. It adds the source term delta*Ab_in and the drug
%                settles at Ab* = delta*Ab_in/(gamma+delta) rather than decaying
%                away. Ab_in = 0 restores the published single pulse.
%   gamma        antibiotic decay rate         (default: 0.0032 for ampicillin,
%                                               0.0050 for colistin)
%   eta          conjugation rate              used when conjugation is on
%   init_As init_Ar init_Bs init_Br init_Cs init_Cr init_P init_L init_M
%                per-compartment initial levels in 1e7 CFU/mL, each overriding
%                whatever the preset and the data would have supplied
%
% So `drug = 0` zeroes epsilon, Ab_init and Ab_in whatever the community says, and
% `drug = 1` restores them from P -- setting epsilon under community 0 no longer
% silently does nothing.
%
% The equations are dYdt_free.m: the published dYdt with the community masking
% dropped, which is what makes the compartments independently settable. The
% masking is redundant for all five presets and this is checked -- see
% tricomm_check_presets.m.
%
% Returned struct:
%   S.t              output times (column)
%   S.Y              nt x 10 raw state [As Ar Bs Br Cs Cr P L M Ab]
%   S.A, S.B, S.C    totals (sensitive + resistant), 1e7 CFU/mL units
%   S.params         the 14-element vector handed to dYdt_free
%   S.y0             initial state actually used
%   S.community, S.antibiotic, S.dose, S.solver
%   S.drug, S.conjugation, S.epsilon, S.eta, S.Ab_init, S.Ab_in, S.gamma
%                    the settings actually in force, after every override
%   S.Ab_star        the plateau the antibiotic settles at, 0 without a feed
%
% ode45 is the default because it is what the published scripts use, so outputs
% reproduce data/model/nodrug_1100h.csv exactly. lsode is ~5x faster and agrees
% to ~1e-4 on these trajectories, which is why the sweeps use it.

if nargin < 2, O = struct(); end

community  = pick(O, P, 'community',  0);
antibiotic = pick(O, P, 'antibiotic', 1);
dose       = pick(O, P, 'dose',       2);
tstart     = pick(O, P, 'tstart',     1);
tmax       = pick(O, P, 'tmax',       NaN);
dt         = pick(O, P, 'dt',         1);
solver     = pick(O, P, 'solver',     'ode45');
data_dir   = pick(O, P, 'data_dir',   fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'model'));

if ~any(community == 0:4)
    error('tricomm_simulate: community must be 0..4, got %g.', community);
end

%% ---- Initial conditions -------------------------------------------------
% Either taken from the experimental time courses (first measured point of the
% selected community, as in tricomm_lin.m) or supplied directly via O.y0.
if isfield(O, 'y0') && ~isempty(O.y0)
    abc = O.y0(:).';
    if numel(abc) ~= 3
        error('tricomm_simulate: y0 must be [A B C], got %d value(s).', numel(abc));
    end
    tdata = [];
else
    [abc, tdata] = ic_from_data(data_dir, community, antibiotic, dose);
end

if ~isfinite(tmax)
    if isempty(tdata)
        error(['tricomm_simulate: tmax is not set and cannot be taken from the ' ...
               'data because y0 was given explicitly. Pass tmax.']);
    end
    tmax = tdata(end);
end
if tmax <= tstart
    error('tricomm_simulate: tmax (%g) must exceed tstart (%g).', tmax, tstart);
end

%% ---- Switches -----------------------------------------------------------
% The community preset, exactly as in the published scripts: the no-drug
% community has no antibiotic and no conjugation; communities 1 and 4 have the
% drug but no conjugation. Explicit drug / conjugation settings win over it.
drug_on = pick(O, P, 'drug',        double(community ~= 0)) ~= 0;
conj_on = pick(O, P, 'conjugation', double(any(community == [2 3]))) ~= 0;

% P.epsilon / P.eta / P.Ab_init are the values used when the switch is on, so
% setting them is meaningful under any community.
if drug_on
    epsilon = pick(O, P, 'epsilon', 1);
    Ab_init = pick(O, P, 'Ab_init', 0.18);
    % Ab_in is the feed concentration, and it FOLLOWS Ab_init unless set on its
    % own. The experiments keep the antibiotic topped up, so one dose number is
    % the honest interface: --ab0 0.32 means "0.32 in the feed, 0.32 at t = 0".
    % --ab-in 0 is the way back to the published single pulse.
    Ab_in   = pick(O, P, 'Ab_in', Ab_init);
else
    epsilon = 0;  Ab_init = 0;  Ab_in = 0;
end
eta = 0;
if conj_on
    eta = pick(O, P, 'eta', 0.005);
end

DEFAULT_GAMMA = [0.0032, 0.0050];    % ampicillin, colistin
gamma = pick(O, P, 'gamma', DEFAULT_GAMMA(antibiotic));

params = [P.mu_A; P.mu_A * P.scale_factor; ...
          P.mu_B; P.mu_B * P.scale_factor; ...
          P.mu_C; P.mu_C * P.scale_factor; ...
          P.alpha_A; P.alpha_B; P.alpha_C; P.delta; epsilon; eta; gamma; Ab_in];

%% ---- Initial state ------------------------------------------------------
% The preset decides whether LM is seeded sensitive or resistant; B and C always
% start fully sensitive, and the amino acids at zero.
switch community
    case {0, 1}
        As0 = abc(1);  Ar0 = 0;
    otherwise
        As0 = 0;       Ar0 = abc(1);
end
y0 = [As0; Ar0; abc(2); 0; abc(3); 0; 0; 0; 0; Ab_init];

% Per-compartment overrides, applied last so they beat both the preset and the
% data. Ab is set from Ab_init above and stays under the drug switch.
INITKEYS = {'init_As','init_Ar','init_Bs','init_Br','init_Cs','init_Cr', ...
            'init_P','init_L','init_M'};
for i = 1:numel(INITKEYS)
    v = pick(O, P, INITKEYS{i}, NaN);
    if isfinite(v)
        y0(i) = v;
    end
end

%% ---- Solve --------------------------------------------------------------
tvec = (tstart:dt:tmax).';
if tvec(end) ~= tmax, tvec(end+1) = tmax; end

f = @(t, Y) dYdt_free(t, Y, params);

% Solver tolerances. Unset means "whatever the solver's own default is", which
% for ode45 is what the published MATLAB scripts ran at and is required to
% reproduce data/model/nodrug_1100h.csv exactly. Set them when the question
% turns on small populations: ode45's default AbsTol of 1e-6 is 10 CFU/mL in
% these units, so extinction studies run at or below the noise floor.
rtol = pick(O, P, 'rtol', NaN);
atol = pick(O, P, 'atol', NaN);

% Absorbing extinction. Without it a strain that falls to 1e-12 CFU/mL is still
% a population the solver can grow back, so a community that really died can
% reappear decades later. extinct_below is in CFU/mL; 0 disables the rule, which
% is the default and is what reproduces the published runs.
extinct_below = pick(O, P, 'extinct_below', 0);
if extinct_below < 0
    error('tricomm_simulate: extinct_below must be >= 0, got %g.', extinct_below);
end
thresh = extinct_below / 1e7;      % CFU/mL -> model units

if thresh > 0
    [t, Y, S_ext] = integrate_absorbing(f, tvec, y0, solver, rtol, atol, thresh);
else
    [t, Y] = integrate_plain(f, tvec, y0, solver, rtol, atol);
    S_ext = struct('strain', {{}}, 'time', []);
end

S = struct();
S.t = t(:);
S.Y = Y;
S.A = Y(:,1) + Y(:,2);
S.B = Y(:,3) + Y(:,4);
S.C = Y(:,5) + Y(:,6);
S.params     = params;
S.y0         = y0;
S.community  = community;
S.antibiotic = antibiotic;
S.dose       = dose;
S.solver     = solver;
S.drug        = drug_on;
S.conjugation = conj_on;
S.epsilon     = epsilon;
S.eta         = eta;
S.Ab_init     = Ab_init;
S.Ab_in       = Ab_in;
S.extinct_below = extinct_below;
S.extinctions   = S_ext;      % .strain (cellstr) and .time, in order
S.gamma       = gamma;
% The plateau the antibiotic relaxes to. Zero without a feed.
S.Ab_star     = Ab_in * P.delta / (gamma + P.delta);
end

% ---------------------------------------------------------------------------
function v = pick(O, P, name, dflt)
% Options struct wins, then the parameter file, then the built-in default.
if isfield(O, name) && ~isempty(O.(name)) && ~(isnumeric(O.(name)) && all(isnan(O.(name))))
    v = O.(name);
elseif isfield(P, name) && ~isempty(P.(name)) && ~(isnumeric(P.(name)) && all(isnan(P.(name))))
    v = P.(name);
else
    v = dflt;
end
end

% ---------------------------------------------------------------------------
function [abc, tdata] = ic_from_data(data_dir, community, antibiotic, dose)
% First measured point of the selected community. Column layout follows
% 1_simulation/linear_scale_version/tricomm_lin.m.
fAMP = fullfile(data_dir, 'pOXA48_interpolated.csv');
fCOL = fullfile(data_dir, 'PN23_interpolated.csv');
if exist(fAMP, 'file') ~= 2
    error(['tricomm_simulate: cannot find %s.\n  Pass data_dir, or pass y0 to ' ...
           'skip the data entirely.'], fAMP);
end

cfu_AMP = readmatrix(fAMP);
if antibiotic == 2
    cfu_COL = readmatrix(fCOL);
    tdata = cfu_COL(:,1);
else
    tdata = cfu_AMP(:,1);
end

switch community
    case 0
        cols = [2 3 4];   src = cfu_AMP;
    case 1
        if antibiotic == 1
            cols = [5 6 7];  src = cfu_AMP;
        else
            src = cfu_COL;
            if dose == 1, cols = [2 3 4]; else, cols = [5 6 7]; end
        end
    case 2
        cols = [8 9 10];
        if antibiotic == 1, src = cfu_AMP; else, src = cfu_COL; end
    case 3
        if antibiotic == 1
            cols = [11 12 13];  src = cfu_AMP;
        else
            src = cfu_COL;
            if dose == 1, cols = [11 12 13]; else, cols = [14 15 16]; end
        end
    case 4
        cols = [14 15 16];  src = cfu_AMP;
end

abc = src(1, cols);
end

% ---------------------------------------------------------------------------
function [t, Y] = integrate_plain(f, tvec, y0, solver, rtol, atol)
% One solver call over the whole grid. ode45 with a two-element tspan returns
% its own internal steps rather than the two requested points, so a midpoint is
% inserted and dropped again -- integrate_absorbing restarts on short spans.
switch lower(solver)
    case 'ode45'
        pad = numel(tvec) == 2;
        tq  = tvec;
        if pad, tq = [tvec(1); mean(tvec); tvec(2)]; end
        if isfinite(rtol) || isfinite(atol)
            oo = odeset();
            if isfinite(rtol), oo = odeset(oo, 'RelTol', rtol); end
            if isfinite(atol), oo = odeset(oo, 'AbsTol', atol); end
            [t, Y] = ode45(f, tq, y0, oo);
        else
            [t, Y] = ode45(f, tq, y0);
        end
        if pad, t = t([1 3]); Y = Y([1 3], :); end
    case 'lsode'
        if ~isfinite(rtol), rtol = 1e-8;  end
        if ~isfinite(atol), atol = 1e-11; end
        lsode_options('relative tolerance', rtol);
        lsode_options('absolute tolerance', atol);
        t = tvec(:);
        Y = lsode(@(Yv, tv) f(tv, Yv), y0, tvec);
    otherwise
        error('tricomm_simulate: unknown solver "%s" (use ode45 or lsode).', solver);
end
end

% ---------------------------------------------------------------------------
function [t, Y, ext] = integrate_absorbing(f, tvec, y0, solver, rtol, atol, thresh)
% Integrate with extinction as an absorbing state: as soon as a strain TOTAL
% (sensitive + resistant) falls below thresh, both of its compartments are set
% to zero and the run continues without them.
%
% Setting the state to zero once is NOT enough. Every term that could refill the
% compartment is proportional to it, so in exact arithmetic it would stay zero --
% but solver roundoff leaves it at ~1e-12 instead, and the growth term amplifies
% that back into a population. So the dead compartments are masked out of the
% right-hand side for the rest of the run: the state is zeroed on the way in and
% the derivative on the way out, which holds them at exactly zero and removes
% their amino acid production, which is the point of the rule.
%
% Detection resolves to the output grid, so a strain is zeroed at the first
% output time at which it is already below thresh -- up to one dt late. Shrink
% --dt if the crossing time itself matters.
PAIRS  = {[1 2], [3 4], [5 6]};
LABELS = {'LM', 'PM', 'PL'};

nt = numel(tvec);
t  = tvec(:);
Y  = zeros(nt, numel(y0));
Y(1, :) = y0(:).';

mask = ones(numel(y0), 1);
dead = false(1, 3);
ext  = struct('strain', {{}}, 'time', []);
i0   = 1;

while i0 < nt
    fm = @(tt, Yv) mask .* f(tt, mask .* Yv);
    [~, Yseg] = integrate_plain(fm, tvec(i0:end), mask .* Y(i0, :).', ...
                                solver, rtol, atol);
    Yseg = Yseg .* mask.';

    hit = 0; hitk = 0;
    for p = 1:3
        if dead(p), continue; end
        k = find(sum(Yseg(:, PAIRS{p}), 2) < thresh, 1);
        if ~isempty(k) && (hitk == 0 || k < hitk)
            hitk = k; hit = p;
        end
    end

    if hit == 0
        Y(i0:nt, :) = Yseg;
        return
    end

    Y(i0:i0+hitk-1, :) = Yseg(1:hitk, :);
    mask(PAIRS{hit}) = 0;
    Y(i0+hitk-1, PAIRS{hit}) = 0;
    dead(hit) = true;
    ext.strain{end+1} = LABELS{hit};
    ext.time(end+1)   = tvec(i0 + hitk - 1);
    i0 = i0 + hitk - 1;      % unchanged when hitk == 1; dead() stops a re-trigger
end
end
