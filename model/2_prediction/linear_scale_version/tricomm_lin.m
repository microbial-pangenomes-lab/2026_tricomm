%% No-drug community simulation compared with data from the long experiment
% This file contains the model simulation of the four no-drug communities
% for 1000 hours and the data of the long experiment.

close all
clear

cfu = readmatrix("../../data/pOXA48_interpolated.csv");
cfu_long = readtable("../../data/pOXA48_long.tsv", "FileType", "text", ...
    "Delimiter", "\t");

%% Experimental data of first experiment: time, means and standard deviations
time = cfu(:,1);
cfu_A = cfu(:,2);
cfu_B = cfu(:,3);
cfu_C = cfu(:,4);

std_A = cfu(:,17);
std_B = cfu(:,18);
std_C = cfu(:,19);

lower_A = cfu(:,20);
lower_B = cfu(:,21);
lower_C = cfu(:,22);

%% Experimental data for long experiment:
cfu_long.time = str2double(erase(string(cfu_long.time), 'h')); % serve per eliminare l'h nei dati
plate_KAN = cfu_long(cfu_long.plate == "KAN", :);
M0_pc = plate_KAN(plate_KAN.vial == "M0", :);
M0_PL_pc = M0_pc(M0_pc.strain == "PL", :);
M0_PM_pc = M0_pc(M0_pc.strain == "PM", :);
M0_LM_pc = M0_pc(M0_pc.strain == "LM", :);

%% Parameters for no-drug communities
mu_A = 0.0003;
mu_B = 0.0005;
mu_C = 0.0004;
alpha_A = 0.7645;
alpha_B = 0.0567;
alpha_C = 0.1202;
delta = 0.0692;
params = [mu_A; mu_B; mu_C; alpha_A; alpha_B; alpha_C; delta];

%% ode45 simulations
% Initial conditions
A_init = cfu_A(1);
B_init = cfu_B(1);
C_init = cfu_C(1);
P_init = 0;
L_init = 0;
M_init = 0;
Y0 = [A_init; B_init; C_init; P_init; L_init; M_init];

[t, Y] = ode45(@(t,Y) dYdt(t,Y,params), 1:1000, Y0);

% Result for each strain and amino acid
A = Y(:, 1);
B = Y(:, 2);
C = Y(:, 3);
P = Y(:, 4);
L = Y(:, 5);
M = Y(:, 6);

% Plot results
figure;
tl = tiledlayout(1,2);

nexttile;
plot(t, C*1e7, 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'PL');
hold on
plot(t, B*1e7, 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'PM');
plot(t, A*1e7, 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.8, 'DisplayName', 'LM');
hold off
xlabel('time (h)', 'FontSize', 15);
ylabel('CFU/mL', 'FontSize', 15);
legend('FontSize', 15,'Box','off');
t1 = title({'Model simulation of bacterial growth for no-drug communities', ''}, ...
    'FontWeight', 'normal');
ax1 = gca;  % get current axes
ax1.FontSize = 16;
ax1.TitleFontSizeMultiplier = 1;   % Deactivate automatic scaling
t1.FontSize = 19;
ax1.XColor = 'black';
ax1.YColor = 'black';
box off

nexttile;
plot(M0_PL_pc.time, M0_PL_pc.CFU_mLAlt, '.--', 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.6, ...
'MarkerSize', 20, 'DisplayName', 'Exp. PL'); hold on
plot(M0_PL_pc.time, M0_PM_pc.CFU_mLAlt, '.--', 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.6, ...
'MarkerSize', 20, 'DisplayName', 'Exp. PM');
plot(M0_PL_pc.time, M0_LM_pc.CFU_mLAlt, '.--', 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.6, ...
'MarkerSize', 20, 'DisplayName', 'Exp. LM');
legend('FontSize', 15,'Box','off');
t2 = title({'Experimental data of long KAN no drug experiment', ''},  'FontWeight', 'normal');
xlabel('time (h)', 'FontSize', 15);
ylabel('CFU/mL', 'FontSize', 15);
ax2 = gca;
ax2.FontSize = 16;
ax2.TitleFontSizeMultiplier = 1; 
t2.FontSize = 19;
ax2.Box = 'off';

linkaxes([ax1 ax2], 'y'); 
ylim(ax1, [2e5 2e8]); 