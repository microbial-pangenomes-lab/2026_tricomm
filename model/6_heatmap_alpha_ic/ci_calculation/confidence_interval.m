%% 95% Confidence Interval Calculation
% In this file, 95% confidence intervals are calculated from the results of 
% the model parameter optimization.

clear; 
clc; 
close all

% Here, CFU non-mLAlt are used, as the SIMPLEXL function needs these values 
% to find the correct parameters. If CFU_mLAlt are used instead, the 
% optimization function does not find the same parameters and the resulting
% fit is different.
cfu = readmatrix("../../data/pOXA48_nonmLAlt_interpolated.csv");

time_A = cfu(:,1);
cfu_A  = cfu(:,2);
cfu_B  = cfu(:,3);
cfu_C  = cfu(:,4);

Data   = [time_A, abs(cfu_A), abs(cfu_B), abs(cfu_C)];

%% SIMPLEXL optimization

p0 = [0.0004, 0.0005, 0.0005, 0.7759, 0.0520, 0.1324, 0.0702];

options = foptions;
options(1)=0; 
options(2)=1e-8; 
options(3)=1e-8; 
options(14)=1000;

[parCal, ~] = SIMPLEXL('objectiveFunction', p0, options, [], Data);
parCal  = parCal(:);
nParams = numel(parCal);

%% Residual vector from objective Function
[~, r0] = objectiveFunction(parCal, Data);
n   = numel(r0);

%% Jacobian matrix calculation
J = zeros(n, nParams);      
h = 1e-6;       
for j = 1:nParams       
    dp = zeros(nParams,1);      
    step = h * max(abs(parCal(j)), 1e-8);   
    dp(j) = step; 

    [~, rPlus]  = objectiveFunction(parCal + dp, Data);
    [~, rMinus] = objectiveFunction(parCal - dp, Data);
    J(:,j) = (rPlus - rMinus) / (2*step);
end

%% 95% Confidence Interval
% The nlparci MATLAB function can be used to calculate the 95% confidence 
% interval.

ci = nlparci(parCal, r0, 'Jacobian', J);
CI_lower = ci(:,1);
CI_upper = ci(:,2);

paramNames = {'mu_A','mu_B','mu_C','alpha_A','alpha_B','alpha_C','delta'};
fprintf('\n%-10s %12s %14s %14s\n', 'Parametro','Stima','IC 95% low','IC 95% up');
for i = 1:nParams
    fprintf('%-10s %12.6g %14.6g %14.6g\n', ...
        paramNames{i}, parCal(i), CI_lower(i), CI_upper(i));
end