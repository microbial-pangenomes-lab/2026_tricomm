function [MSE, residual] = objectiveFunction(Parameters,Data)

vTime = Data(:,1);
cfu_A = Data(:,2);
cfu_B = Data(:,3);
cfu_C = Data(:,4);

% Initial points
A_init = cfu_A(1);
B_init = cfu_B(1);
C_init = cfu_C(1);
P_init = 0;
L_init = 0;
M_init = 0;

y0 = [A_init; B_init; C_init; P_init; L_init; M_init];

mu_A = Parameters(1);
mu_B = Parameters(2);
mu_C = Parameters(3);
alpha_A = Parameters(4);
alpha_B = Parameters(5);
alpha_C = Parameters(6);
delta   = Parameters(7);

params = [mu_A; mu_B; mu_C; alpha_A; alpha_B; alpha_C; delta];

% ODE
[tsol,sol] = ode45(@(t,y) dYdt(t,y,params), vTime, y0); % da inserire la function (modello da calibrare)

% Results
solCfuA=sol(:,1);
solCfuB=sol(:,2);
solCfuC=sol(:,3);

% Plot creation 
% figure(1)
% plot(tsol,solCfuA,'m', vTime, cfu_A, '*m', ...
%      tsol,solCfuB,'r', vTime, cfu_B, '*r', ...
%      tsol,solCfuC,'g', vTime, cfu_C, '*g')
% drawnow

% Method for calculating the error
method=1;

% Squared error
errA = (solCfuA - cfu_A).^2;
errB = (solCfuB - cfu_B).^2;
errC = (solCfuC - cfu_C).^2;

switch method
    case 1
        MSE=(1/length(vTime))*sum(errA+errB+errC);
    case 2
        MSE=(1/length(vTime))*sum((cfu_A.*(cfu_A-solCfuA)).^2 + ...
            (cfu_B.*(cfu_B-solCfuB)).^2 + (cfu_C.*(cfu_C-solCfuC)).^2);

end

residual = [solCfuA - cfu_A; solCfuB - cfu_B; solCfuC - cfu_C];


