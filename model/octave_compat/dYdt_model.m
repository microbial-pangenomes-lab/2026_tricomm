function Ydot = dYdt_model(t, Y, params, community, antibiotic)
% DYDT_MODEL  10-state tricomm ODE (drug / no-drug communities).
%   Identical to 1_simulation/linear_scale_version/dYdt.m, renamed so it
%   never shadows the repo's own dYdt.m files when octave_compat is on the
%   Octave path. Used by run_model.m.
As = Y(1);
Ar = Y(2);
Bs = Y(3);
Br = Y(4);
Cs = Y(5);
Cr = Y(6);
P = Y(7);
L = Y(8);
M = Y(9);
Ab = Y(10);

mu_A = params(1);
mu_Ar = params(2);
mu_B = params(3);
mu_Br = params(4);
mu_C = params(5);
mu_Cr = params(6);
alpha_A = params(7);
alpha_B = params(8);
alpha_C = params(9);
delta = params(10);
epsilon = params(11);
eta = params(12);

switch antibiotic
    case 1
        gamma = 0.0032;
    case 2
        gamma = 0.0050;
end

switch community
    case {0,1}
        Ar = 0;
        Br = 0;
        Cr = 0;
    case {2,3}
        As = 0;
    case 4
        As = 0;
        Br = 0;
        Cr = 0;
end

dAs_dt = mu_A * As * L * M - delta * As - epsilon * Ab * As;
dAr_dt = mu_Ar * Ar * L * M - delta * Ar;
dBs_dt = mu_B * Bs * P * M - delta * Bs - epsilon * Ab * Bs - eta * Bs * (Ar + Br + Cr);
dBr_dt = mu_Br * Br * P * M - delta * Br + eta * Bs * (Ar + Br + Cr);
dCs_dt = mu_C * Cs * P * L - delta * Cs - epsilon * Ab * Cs - eta * Cs * (Ar + Br + Cr);
dCr_dt = mu_Cr * Cr * P * L - delta * Cr + eta * Cs * (Ar + Br + Cr);
dP_dt = -mu_B * (Bs + Br) * P * M - mu_C * (Cs + Cr) * P * L + alpha_A * (As + Ar);
dL_dt = -mu_A * (As + Ar) * L * M - mu_C * (Cs + Cr) * P * L + alpha_B * (Bs + Br);
dM_dt = -mu_A * (As + Ar) * L * M - mu_B * (Bs + Br) * P * M + alpha_C * (Cs + Cr);
dAb_dt = -(gamma + delta) * Ab;

Ydot = [dAs_dt; dAr_dt; dBs_dt; dBr_dt; dCs_dt; dCr_dt; dP_dt; dL_dt; dM_dt; dAb_dt];
end
