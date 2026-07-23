function Ydot = dYdt(t, Y, params)
Ar = Y(1);
Bs = Y(2);
Br = Y(3);
Cs = Y(4);
Cr = Y(5);
P = Y(6);
L = Y(7);
M = Y(8);
Ab = Y(9);

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

gamma = 0.0032;

dAr_dt = mu_Ar * Ar * L * M - delta * Ar;
dBs_dt = mu_B * Bs * P * M - delta * Bs - epsilon * Ab * Bs - eta * Bs * (Ar + Br + Cr);
dBr_dt = mu_Br * Br * P * M - delta * Br + eta * Bs * (Ar + Br + Cr);
dCs_dt = mu_C * Cs * P * L - delta * Cs - epsilon * Ab * Cs - eta * Cs * (Ar + Br + Cr);
dCr_dt = mu_Cr * Cr * P * L - delta * Cr + eta * Cs * (Ar + Br + Cr);
dP_dt = -mu_B * (Bs + Br) * P * M - mu_C * (Cs + Cr) * P * L + alpha_A * Ar;
dL_dt = -mu_A * Ar * L * M - mu_C * (Cs + Cr) * P * L + alpha_B * (Bs + Br);
dM_dt = -mu_A * Ar * L * M - mu_B * (Bs + Br) * P * M + alpha_C * (Cs + Cr);
dAb_dt = -(gamma + delta) * Ab;

Ydot = [dAr_dt; dBs_dt; dBr_dt; dCs_dt; dCr_dt; dP_dt; dL_dt; dM_dt; dAb_dt];
end