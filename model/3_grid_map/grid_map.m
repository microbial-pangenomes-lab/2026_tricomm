%% Grid map realization for third drug community
% In this file a grid map is created with plots of different values of 
% epsilon and eta. 

close all
clear
clc

cfu = readmatrix("../data/pOXA48_interpolated.csv");

% n indicated the number of eta and epsilon that are considered.
n = 10;
% linspace creates a series of values as linspace(starting_value, 
% final_value, number_of_desired_values)
epsilon_values = linspace(0.1,1,n);
eta_values = linspace(0.00001,0.001,n);

% scale_factor for mu_A, mu_B and mu_C
scale_factor = 1.1;

ax_sum = cell(n,n);   % axes for the sum of strains
ax_sep = cell(n,n);   % axes for the separated strains

% Possible which_plot:
% 'sum' = sum of resistant and sensitive strains
% 'separated' = separated resistant and sensitive strains

which_plot = 'separated';

for i = 1:n
    epsilon = epsilon_values(i);
    for k = 1:n
        eta = eta_values(k);
        
        %% Choice of experimental data
        % In this part, experimental data are selected with respect to which
        % community is observed.

        time = cfu(:,1);

        cfu_A = cfu(:, 11);
        cfu_B = cfu(:, 12);
        cfu_C = cfu(:, 13);
        
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
        
        params = [mu_A; mu_Ar; mu_B; mu_Br; mu_C; mu_Cr; alpha_A; alpha_B; alpha_C; 
            delta; epsilon; eta];
        
        %% ode45 simulations
        % Initial conditions
        Ar_init = cfu_A(1);
        Bs_init = cfu_B(1);
        Br_init = 0;
        Cs_init = cfu_C(1);
        Cr_init = 0;
        P_init = 0;
        L_init = 0;
        M_init = 0;
        Ab_init = 0.18;
        Y0 = [Ar_init; Bs_init; Br_init; Cs_init; Cr_init; P_init; L_init; 
            M_init; Ab_init];
        
        [t, Y] = ode45(@(t,Y) dYdt(t,Y,params), 1:time(end), Y0);
        
        % Result for each strain and amino acid
        Ar = Y(:, 1);
        Bs = Y(:, 2);
        Br = Y(:, 3);
        Cs = Y(:, 4);
        Cr = Y(:, 5);
        P = Y(:, 6);
        L = Y(:, 7);
        M = Y(:, 8);
        Ab = Y(:, 9);
        
        B = Bs + Br;
        C = Cs + Cr;

        %% Plot with the sum of sensitive and resistant strains
        figure('Visible','off');
        plot(t, C, 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'PL'); hold on
        plot(t, B, 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'PM');
        plot(t, Ar, 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.8, 'DisplayName', 'LM'); 
        xlabel('time (h)', 'FontSize', 15);
        ylabel('CFU/mL', 'FontSize', 15);
        legend('FontSize', 15, 'Box','off');
        sgtitle('Model simulation of bacterial growth for PVI 16xMIC', ...
            'FontSize', 19);
        ax.FontSize = 16;
        ax.XColor = 'black';
        ax.YColor = 'black';
        ax_sum{i,k} = ax;
        
        %% Plot with separated sensitive and resistant strains
        figure('Visible', 'off');
        plot(t, Cs, 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'Sensitive PL'); hold on
        plot(t, Cr, '--', 'Color', [0.20, 0.75, 0.60], 'LineWidth', 1.8, 'DisplayName', 'Resistant PL');
        plot(t, Bs, 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'Sensitive PM');
        plot(t, Br, '--', 'Color', [0.85, 0.25, 0.25], 'LineWidth', 1.8, 'DisplayName', 'Resistant PM');
        plot(t, Ar, 'Color', [0.75, 0.75, 0.75], 'LineWidth', 1.8, 'DisplayName', 'Resistant LM'); 
        xlabel('time (h)', 'FontSize', 15);
        ylabel('CFU/mL', 'FontSize', 15);
        legend('FontSize', 15, 'Box','off');
        sgtitle('Model simulation of bacterial growth for PVI 16xMIC', ...
                'FontSize', 19)
        ax = gca;  % get current axes
        ax.FontSize = 16;
        ax.XColor = 'black';
        ax.YColor = 'black';
        box off
        ax_sep{i,k} = gca;

    end

end


%% Plots
switch which_plot
    case 'sum'
        source_axes = ax_sum;
    case 'separated'
        source_axes = ax_sep;
end

grid_fig = figure;
tl = tiledlayout(grid_fig, n, n, 'TileSpacing', 'compact', 'Padding', 'compact'); % This contains 25 void cells

for i = 1:n
    for k = 1:n
        src_ax = source_axes{i, k};

        nexttile(tl);       % This makes the new cell active -> every call is in a new cell
        dest_ax = gca;
        copyobj(allchild(src_ax), dest_ax);   % This part copies all the part inside src_ax and put them in dest_ax
        dest_ax.XLim = src_ax.XLim;
        dest_ax.YLim = src_ax.YLim;

        if k == 1
            ylabel(dest_ax, sprintf('\\epsilon = %.4g', epsilon_values(i)), 'FontSize', 12, ...
                'Rotation', 0);
        end
        if i == 1
            title(dest_ax, sprintf('\\eta = %.4g', eta_values(k)), 'FontSize', 12, ...
                'FontWeight', 'normal');
        end
        if i == n
            xlabel(dest_ax, 'time (h)', 'FontSize', 12);
        end
    end
end
lgd = legend(dest_ax);
lgd.Layout.Tile = 'east';
lgd.Orientation = 'vertical';
lgd.FontSize = 14;
lgd.Box = 'off';
title(tl, 'Grid map for \epsilon and \eta parameters', 'FontSize', 18);
