function Ydot = dYdt_free(t, Y, params)
% DYDT_FREE  The 10-state TriComm ODE with nothing forced to zero.
%
%   Ydot = dYdt_free(t, Y, params)
%   params = [mu_A mu_Ar mu_B mu_Br mu_C mu_Cr alpha_A alpha_B alpha_C ...
%             delta epsilon eta gamma Ab_in]
%
% Same equations as dYdt_model.m (the unmodified published dYdt), with three
% differences:
%
%   * no `community` switch. dYdt_model zeroes compartments per community --
%     Ar/Br/Cr for communities 0 and 1, As for 2 and 3, As/Br/Cr for 4. That
%     masking is redundant: every term of dAs, dAr, dBr, dCr that could seed a
%     compartment from zero is either proportional to that compartment or
%     carries the factor eta, so a compartment that starts at zero with the
%     matching rate off stays at zero on its own. Dropping the switch is what
%     lets the caller set the compartments independently of the preset, and
%     tricomm_simulate verifies the two agree on all five presets.
%
%   * gamma (antibiotic decay) comes in through params rather than an
%     `antibiotic` code, so the dose regime is a number the caller chooses.
%
%   * Ab_in, the antibiotic concentration in the feed, adds the source term
%     delta * Ab_in. The published model has no source: the drug is a single
%     pulse decaying at gamma + delta, so it is effectively gone after ~3
%     dilution times whatever the experiment does. With a feed the drug relaxes
%     to a plateau instead,
%
%         Ab(t) -> Ab* = delta * Ab_in / (gamma + delta)
%
%     approached with the same time constant 1/(gamma + delta). Ab_in = 0
%     reproduces the published pulse exactly.
%
%     The callers default Ab_in to Ab_init, because that is what the experiments
%     do: the drug is held at the dose rather than allowed to wash out. The run
%     then starts at the feed concentration and settles just below it, short by
%     exactly the fraction gamma degrades, gamma / (gamma + delta).
%
% Units: strains in 1e7 CFU/mL, time in hours.

As = Y(1);  Ar = Y(2);
Bs = Y(3);  Br = Y(4);
Cs = Y(5);  Cr = Y(6);
P  = Y(7);  L  = Y(8);  M = Y(9);  Ab = Y(10);

mu_A = params(1);  mu_Ar = params(2);
mu_B = params(3);  mu_Br = params(4);
mu_C = params(5);  mu_Cr = params(6);
alpha_A = params(7);  alpha_B = params(8);  alpha_C = params(9);
delta   = params(10); epsilon = params(11); eta     = params(12);
gamma   = params(13);
% Optional 14th entry so a 13-element params vector still runs as the published
% pulse rather than erroring.
if numel(params) >= 14, Ab_in = params(14); else, Ab_in = 0; end

% Conjugation: any resistant cell can hand the plasmid to a sensitive B or C.
donors = Ar + Br + Cr;

dAs_dt = mu_A  * As * L * M - delta * As - epsilon * Ab * As;
dAr_dt = mu_Ar * Ar * L * M - delta * Ar;
dBs_dt = mu_B  * Bs * P * M - delta * Bs - epsilon * Ab * Bs - eta * Bs * donors;
dBr_dt = mu_Br * Br * P * M - delta * Br                     + eta * Bs * donors;
dCs_dt = mu_C  * Cs * P * L - delta * Cs - epsilon * Ab * Cs - eta * Cs * donors;
dCr_dt = mu_Cr * Cr * P * L - delta * Cr                     + eta * Cs * donors;

dP_dt = -mu_B * (Bs + Br) * P * M - mu_C * (Cs + Cr) * P * L + alpha_A * (As + Ar);
dL_dt = -mu_A * (As + Ar) * L * M - mu_C * (Cs + Cr) * P * L + alpha_B * (Bs + Br);
dM_dt = -mu_A * (As + Ar) * L * M - mu_B * (Bs + Br) * P * M + alpha_C * (Cs + Cr);
dAb_dt = delta * Ab_in - (gamma + delta) * Ab;

Ydot = [dAs_dt; dAr_dt; dBs_dt; dBr_dt; dCs_dt; dCr_dt; dP_dt; dL_dt; dM_dt; dAb_dt];
end
