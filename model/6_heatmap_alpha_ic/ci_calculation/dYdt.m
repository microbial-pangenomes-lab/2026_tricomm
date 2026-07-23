function Ydot = dYdt(t, Y, params)

A = Y(1);
B = Y(2);
C = Y(3);
P = Y(4);
L = Y(5);
M = Y(6);

mu_A = params(1);
mu_B = params(2);
mu_C = params(3);
alpha_A = params(4);
alpha_B = params(5);
alpha_C = params(6);
delta = params(7);

dA_dt = mu_A * A * L * M - delta * A;
dB_dt = mu_B * B * P * M - delta * B;
dC_dt = mu_C * C * P * L - delta * C;
dP_dt = -mu_B * B * P * M - mu_C * C * P * L + alpha_A * A;
dL_dt = -mu_A * A * L * M - mu_C * C * P * L + alpha_B * B;
dM_dt = -mu_A * A * L * M - mu_B * B * P * M + alpha_C * C;
Ydot = [dA_dt; dB_dt; dC_dt; dP_dt; dL_dt; dM_dt];

end