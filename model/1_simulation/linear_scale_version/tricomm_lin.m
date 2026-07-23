%% Final model for Tricomm
% This file contains the model of the drug and no-drug communities and 
% uses the parameters derived from the no-drug communities. 

close all
clear
clc

cfu_AMP = readmatrix("../../data/pOXA48_interpolated.csv");
cfu_COL = readmatrix("../../data/PN23_interpolated.csv");

%% Choise of community
% It is possible to choose which antibiotic experiments to visualize, indicating
% in "antibiotic" variable:
% 1 = ampicillin 16xMIC
% 2 = colistin
% For colistin, you can choose the "dose" of the communities 1 and 3 between
% 1 = 1xMIC
% 2 = 2xMIC

antibiotic = 1;

if antibiotic == 2
    dose = 2;
end

% It is possible to choose also  which community to visualize, indicating the 
% number in "community" variable:
% For AMPICILLIN EXPERIMENTS
% 0 = no-drug community
% 1 = first drug community -> all strains are sensitive
% 2 = second drug community -> all strains are sensitive and after two
% hours a resistent LM strains is added to the community. It is able to 
% transfer its plasmid
% 3 = third drug community -> LM is resistent and it transfers its plasmid
% 4 = fourth drug community -> LM is resistent and it does not transfer its
% plasmid

% For COLISTIN EXPERIMENT
% 0 = no-drug community -> same values of ampicillin experiment
% 1 = KAN drug
% 2 = KAN + PN23
% 3 = PN23 drug

community = 3;

% All community parameters are automatically modified just by changing the
% community number. epsilon, eta, Ab_init and scale_factor can be chosen 
% here. epsilon is the antibiotic-induced killing rate, eta is the
% conjugation rate for plasmid transfer, Ab_init is the initial
% concentration of antibiotic and scale_factor is a scaling factor for the
% growth rate of resistant strains, considering 1 as the growth rate of
% sensitive strains. It includes plasmid costs to bacterial fitness.

epsilon = 1;
eta = 0.005;
Ab_init = 0.18;
scale_factor = 1.1;

%% Choice of experimental data
% In this part, experimental data are selected with respect to which
% community is observed.

switch antibiotic
    case 1
        time = cfu_AMP(:,1);
    case 2
        time = cfu_COL(:,1);
end

switch community
            case 0      % no-drug community
                cfu_A = cfu_AMP(:, 2);
                cfu_B = cfu_AMP(:, 3);
                cfu_C = cfu_AMP(:, 4);

            case 1      % first community
                switch antibiotic
                    case 1
                        cfu_A = cfu_AMP(:, 5);
                        cfu_B = cfu_AMP(:, 6);
                        cfu_C = cfu_AMP(:, 7);
                    case 2
                        switch dose
                            case 1
                                cfu_A = cfu_COL(:, 2);
                                cfu_B = cfu_COL(:, 3);
                                cfu_C = cfu_COL(:, 4);
                            case 2
                                cfu_A = cfu_COL(:, 5);
                                cfu_B = cfu_COL(:, 6);
                                cfu_C = cfu_COL(:, 7);
                        end
                end
            case 2      % second community
                switch antibiotic
                    case 1
                        cfu_A = cfu_AMP(:, 8);
                        cfu_B = cfu_AMP(:, 9);
                        cfu_C = cfu_AMP(:, 10);
                    case 2
                        cfu_A = cfu_COL(:, 8);
                        cfu_B = cfu_COL(:, 9);
                        cfu_C = cfu_COL(:, 10);
                end
            case 3      % third community
                switch antibiotic
                    case 1
                        cfu_A = cfu_AMP(:, 11);
                        cfu_B = cfu_AMP(:, 12);
                        cfu_C = cfu_AMP(:, 13);
                    case 2
                        switch dose
                            case 1
                                cfu_A = cfu_COL(:, 11);
                                cfu_B = cfu_COL(:, 12);
                                cfu_C = cfu_COL(:, 13);
                            case 2
                                cfu_A = cfu_COL(:, 14);
                                cfu_B = cfu_COL(:, 15);
                                cfu_C = cfu_COL(:, 16);
                        end
                end
            case 4      % fourth community
                cfu_A = cfu_AMP(:, 14);
                cfu_B = cfu_AMP(:, 15);
                cfu_C = cfu_AMP(:, 16);
end

% Standard deviation of no-drug community from ampicillin experiments
std_A = cfu_AMP(:,17);
std_B = cfu_AMP(:,18);
std_C = cfu_AMP(:,19);

% Lower limit for the standard deviation of the no-drug communities from 
% ampicillin experiments. This is important to set the lower error bar to 
% zero if it is negative.
lower_A = cfu_AMP(:,20);
lower_B = cfu_AMP(:,21);
lower_C = cfu_AMP(:,22);

%% Parameters from no-drug communities of ampicillin experiments
mu_A = 0.0003;
mu_Ar = mu_A * scale_factor;
mu_B = 0.0005;
mu_Br = mu_B * scale_factor;
mu_C = 0.0004;
mu_Cr = mu_C * scale_factor;
alpha_A = 0.7645;
alpha_B = 0.0567;
alpha_C = 0.1202;
delta = 0.0692;

if community == 0
    epsilon = 0;
    eta = 0;
    Ab_init = 0;
end

if community == 1 || community == 4
    eta = 0;
end

params = [mu_A; mu_Ar; mu_B; mu_Br; mu_C; mu_Cr; alpha_A; alpha_B; alpha_C; 
    delta; epsilon; eta];

%% ode45 simulations
% Initial conditions
switch community
    case {0, 1}
        As_init = cfu_A(1);
        Ar_init = 0;
    case {2, 3, 4}
        As_init = 0;
        Ar_init = cfu_A(1);
end
Bs_init = cfu_B(1);
Br_init = 0;
Cs_init = cfu_C(1);
Cr_init = 0;
P_init = 0;
L_init = 0;
M_init = 0;
Y0 = [As_init; Ar_init; Bs_init; Br_init; Cs_init; Cr_init; P_init; L_init; 
    M_init; Ab_init];

[t, Y] = ode45(@(t,Y) dYdt(t,Y,params, community, antibiotic), 1:time(end), Y0);

% Result for each strain and amino acid
As = Y(:, 1);
Ar = Y(:, 2);
Bs = Y(:, 3);
Br = Y(:, 4);
Cs = Y(:, 5);
Cr = Y(:, 6);
P = Y(:, 7);
L = Y(:, 8);
M = Y(:, 9);
Ab = Y(:, 10);

B = Bs + Br;
C = Cs + Cr;
A = As + Ar;

% Community names
if antibiotic == 1
    switch community
        case 0
            community_name = "no-drug community";
        case 1
            community_name = "KAN";
        case 2
            community_name = "KAN + PVI";
        case 3
            community_name = "PVI";
        case 4
            community_name = "PVB";
    end
else
    switch community
        case 0
            community_name = "no-drug community";
        case 1
            community_name = "KAN";
        case 2
            community_name = "KAN + PN23";
        case 3
            community_name = "PN23";
    end
end

switch antibiotic
    case 1
        antibiotic_name = "16xMIC";
    case 2
        if community == 1 || community == 3
            switch dose
                case 1
                    antibiotic_name = "1xMIC";
                case 2
                    antibiotic_name = "2xMIC";
            end
        elseif community == 2
            antibiotic_name = "2xMIC";
        end
end

%% Plot with the sum of sensitive and resistant strains
figure;
tl = tiledlayout(1,2);

if community == 0
    nexttile;
    plot(t, C*1e7, 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'Fit PL'); hold on
    plot(t, B*1e7, 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'Fit PM');
    plot(t, A*1e7, 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.8, 'DisplayName', 'Fit LM'); 
    e1 = errorbar(time, cfu_C*1e7, lower_C*1e7, std_C*1e7, '.', 'Color', [0.20, 0.75, 0.60], 'LineWidth', 0.85, ...
        'MarkerSize', 20, 'DisplayName', 'Exp. PL');
    e2 = errorbar(time, cfu_B*1e7, lower_B*1e7, std_B*1e7, '.', 'Color', [0.85, 0.25, 0.25], 'LineWidth', 0.85, ...
      'MarkerSize', 20, 'DisplayName', 'Exp. PM');
    e3 = errorbar(time, cfu_A*1e7, lower_A*1e7, std_A*1e7, '.', 'Color', [0.75, 0.75, 0.75], 'LineWidth', 0.85, ...
        'MarkerSize', 20, 'DisplayName', 'Exp. LM');
    e1.CapSize = 0;
    e2.CapSize = 0;
    e3.CapSize = 0;
    
    xlabel('time (h)', 'FontSize', 15);
    ylabel('CFU/mL', 'FontSize', 15);
    legend('FontSize', 15, 'Box','off');
    t1 = title({'Model fitting of bacterial growth for no-drug community', ''}, ...
    'FontWeight', 'normal');
    ax = gca;  % get current axes
    ax.FontSize = 16;
    ax.TitleFontSizeMultiplier = 1;   % Deactivate automatic scaling
    t1.FontSize = 19;
    ax.XColor = 'black';
    ax.YColor = 'black';
    box off
else
    % Main plot
    ax_main = nexttile;
    hold(ax_main, 'on')
    
    plot(ax_main, t, C*1e7, 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'PL'); hold on
    plot(ax_main, t, B*1e7, 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'PM');
    plot(ax_main, t, A*1e7, 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.8, 'DisplayName', 'LM'); 
    
    y_max_first200 = max([C(20:200); B(20:200); A(20:200)]) *1e7;   % It is useful not to overlap main plot with little plot
    hold(ax_main, 'off');
    xlabel(ax_main, 'time (h)', 'FontSize', 15);
    ylabel(ax_main, 'CFU/mL', 'FontSize', 15);
    ylim(ax_main, [0, y_max_first200 + 15*1e7]);
    legend(ax_main, 'FontSize', 15, 'Box','off', 'Location', 'southoutside', 'Orientation','horizontal');
    t2 = title({sprintf('Model simulation of bacterial growth for %s %s', community_name, antibiotic_name), ''}, ...
            'FontWeight', 'normal', 'FontSize', 19);
    ax_main.FontSize = 16;
    ax_main.TitleFontSizeMultiplier = 1;  
    t2.FontSize = 19;
    ax_main.XColor = 'black';
    ax_main.YColor = 'black';
    box(ax_main, 'off');
end

%% Plot with separated sensitive and resistant strains

if community ~= 0
    nexttile;
    switch community
        case {1}
            plot(t, Cs*1e7, 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'Sensitive PL'); hold on
            plot(t, Bs*1e7, 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'Sensitive PM');
            plot(t, As*1e7, 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.8, 'DisplayName', 'Sensitive LM'); 
        case {2, 3}
            plot(t, Cs*1e7, 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'Sensitive PL'); hold on
            plot(t, Cr*1e7, '--', 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'Resistant PL');
            plot(t, Bs*1e7, 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'Sensitive PM');
            plot(t, Br*1e7, '--', 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'Resistant PM');
            plot(t, Ar*1e7, 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.8, 'DisplayName', 'Resistant LM'); 
        case 4
            plot(t, Cs*1e7, 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'Sensitive PL'); hold on
            plot(t, Bs*1e7, 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'Sensitive PM');
            plot(t, Ar*1e7, 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.8, 'DisplayName', 'Resistant LM'); 
    end
    xlabel('time (h)', 'FontSize', 15);
    ylabel('CFU/mL', 'FontSize', 15);
    if community ~= 0       % For drug communities
        ylim([0, y_max_first200 + 15*1e7]);
    end
    legend('FontSize', 15, 'Box','off', 'Location', 'southoutside', 'Orientation','horizontal');
    t3 = title({sprintf('Model simulation of bacterial growth for %s %s', community_name, antibiotic_name), ''}, ...
            'FontWeight', 'normal', 'FontSize', 19);
    ax = gca;  % get current axes
    ax.FontSize = 16;
    ax.TitleFontSizeMultiplier = 1;   % Deactivate automatic scaling
    t3.FontSize = 19;
    ax.XColor = 'black';
    ax.YColor = 'black';
    box off
end

if community ~= 0
    % Little plot inside the main one
    pos = ax_main.Position; 
    ax_inset = axes('Position', [pos(1)-0.02, pos(2)+0.45, 0.38*pos(3), 0.3*pos(4)]);
    hold(ax_inset, 'on');
    
    plot(ax_inset, time, cfu_C*1e7, '.', 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.6, ...
    'MarkerSize', 20, 'DisplayName', 'Exp. PL');
    plot(ax_inset, time, cfu_B*1e7, '.', 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.6, ...
    'MarkerSize', 20, 'DisplayName', 'Exp. PM');
    plot(ax_inset, time, cfu_A*1e7, '.', 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.6, ...
    'MarkerSize', 20, 'DisplayName', 'Exp. LM');

    hold(ax_inset, 'off');
    title(ax_inset, 'Experimental data');
    ax_inset.FontSize = 9;
    ax_inset.Title.FontSize = 12;
    ax_inset.Box = 'on';
    ax_inset.YScale = 'log';
end

%--------------------------------------------------------------------%
%% Simulation after 306 hours for no-drug communities
% In this part, the model simulate oscillating trend over 1000 hours. 

if community == 0
    % 1:1000 indicates that the simulation runs up to 1000 hours.
    [t, Y] = ode45(@(t,Y) dYdt(t,Y,params, community, antibiotic), 1:1000, Y0);

    % Result for each strain and amino acid
    As = Y(:, 1);
    Ar = Y(:, 2);
    Bs = Y(:, 3);
    Br = Y(:, 4);
    Cs = Y(:, 5);
    Cr = Y(:, 6);
    P = Y(:, 7);
    L = Y(:, 8);
    M = Y(:, 9);
    Ab = Y(:, 10);

    % Plot results
    nexttile;
    plot(t, Cs*1e7, 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'PL');
    hold on
    plot(t, Bs*1e7, 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'PM');
    plot(t, As*1e7, 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.8, 'DisplayName', 'LM');
    hold off
    xlabel('time (h)', 'FontSize', 15);
    ylabel('CFU/mL', 'FontSize', 15);
    legend('FontSize', 15,'Box','off');
    t1 = title({'Model simulation of bacterial growth for no-drug communities', ''}, ...
        'FontWeight', 'normal');
    ax = gca;  % get current axes
    ax.FontSize = 16;
    ax.TitleFontSizeMultiplier = 1;   % Deactivate automatic scaling
    t1.FontSize = 19;
    ax.XColor = 'black';
    ax.YColor = 'black';
    box off
end