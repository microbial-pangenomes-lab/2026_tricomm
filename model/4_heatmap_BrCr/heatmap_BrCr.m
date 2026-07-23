%% Heatmap realization for third drug community
% In this file a heatmap is created with different values of epsilon and
% eta. 

close all
clear
clc

cfu = readmatrix("../data/pOXA48_interpolated.csv");

% n indicated the number of eta and epsilon that are considered.
n = 25;
% linspace creates a series of values as linspace(starting_value, 
% final_value, number_of_desired_values)
epsilon_values = linspace(0.1,1,n);
eta_values = linspace(0.0005,0.0055,n);

% scale_factor for mu_A, mu_B and mu_C
scale_factor = 1.1;

% Table creation for time values
crossing_time_1 = nan(n,n);
crossing_time_2 = nan(n,n);
crossing_time_3 = nan(n,n);
Br_Btot_time_1  = nan(n,n);
Cr_Ctot_time_1  = nan(n,n);

for i = 1:n
    epsilon = epsilon_values(i);
    for k = 1:n
        eta = eta_values(k);

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
        Br_Btot = Br./B;
        Cr_Ctot = Cr./C;

        % Search for the first time point at which both Br and Cr 
        % exceed 3 CFU, Br is the 50% of the total B strain and Cr is
        % the 50% of the total C strain.
        idx_cross = find(Br > 3 & Cr > 3, 1, 'first');
        idx_cross_Br = find(Br > 3, 1, 'first');
        idx_cross_Cr = find(Cr > 3, 1, 'first');
        idx_Br    = find(Br_Btot > 0.5, 1, 'first');
        idx_Cr    = find(Cr_Ctot > 0.5, 1, 'first');
        
        if ~isempty(idx_cross), crossing_time_1(i,k) = t(idx_cross); end
        if ~isempty(idx_cross_Br), crossing_time_2(i,k) = t(idx_cross_Br); end
        if ~isempty(idx_cross_Cr), crossing_time_3(i,k) = t(idx_cross_Cr); end
        if ~isempty(idx_Br),    Br_Btot_time_1(i,k)  = t(idx_Br);    end
        if ~isempty(idx_Cr),    Cr_Ctot_time_1(i,k)  = t(idx_Cr);    end
    end

end

%% Heatmap plot (epsilon on columns, eta on rows)
epsilon_labels = arrayfun(@(x) sprintf('%.4g', x), epsilon_values, 'UniformOutput', false);
eta_labels = arrayfun(@(x) sprintf('%.4g', x), eta_values, 'UniformOutput', false);

figure;
heatmap_span = 5;   % Grid columns occupied by every heatmap
spacer_span  = 1;   % Void grid column between the heatmaps
ncols = 3*heatmap_span + 2*spacer_span;

tl = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
tl.OuterPosition = [0 0 0.96 1];

tl_row1 = tiledlayout(tl, 1, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');
tl_row1.Layout.Tile = 1;
tl_row2 = tiledlayout(tl, 1, ncols, 'TileSpacing', 'compact', 'Padding', 'compact');
tl_row2.Layout.Tile = 2;

col1 = 1;
col2 = col1 + heatmap_span + spacer_span;
col3 = col2 + heatmap_span + spacer_span;

%% Heatmap of crossing times for Br AND Cr (epsilon on columns, eta on rows)
nexttile(tl_row1, col1, [1 heatmap_span]);
h1 = heatmap(epsilon_labels, eta_labels, crossing_time_1', 'Colormap', parula);
h1.FontSize = 13;
h1.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h1.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h1.Title = 'Br & Cr > 3 CFU';
h1.XLabel = '\epsilon';
h1.YLabel = '\eta';
h1.CellLabelColor = 'none';
h1.MissingDataColor = [1 1 1];
if any(isnan(h1.ColorData(:)))
    h1.MissingDataLabel = 'no resistant';
end

%% Heatmap of crossing times for Br (epsilon on columns, eta on rows)
nexttile(tl_row1, col2, [1 heatmap_span]);
h2 = heatmap(epsilon_labels, eta_labels, crossing_time_2', 'Colormap', parula);
h2.FontSize = 13;
h2.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h2.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h2.Title = 'Br > 3 CFU';
h2.XLabel = '\epsilon';
h2.YLabel = '\eta';
h2.CellLabelColor = 'none';
h2.MissingDataColor = [1 1 1];
if any(isnan(h2.ColorData(:)))
    h2.MissingDataLabel = 'no resistant';
end

%% Heatmap of crossing times for Cr (epsilon on columns, eta on rows)
nexttile(tl_row1, col3, [1 heatmap_span]);
h3 = heatmap(epsilon_labels, eta_labels, crossing_time_3', 'Colormap', parula);
h3.FontSize = 13;
h3.Title = 'Cr > 3 CFU';
h3.XLabel = '\epsilon';
h3.YLabel = '\eta';
h3.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h3.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h3.CellLabelColor = 'none';
h3.MissingDataColor = [1 1 1];
if any(isnan(h3.ColorData(:)))
    h3.MissingDataLabel = 'no resistant';
end

%% Heatmap of crossing times for Br/Btot (epsilon on columns, eta on rows)
nexttile(tl_row2, col1, [1 heatmap_span]);
h4 = heatmap(epsilon_labels, eta_labels, Br_Btot_time_1', 'Colormap', parula);
h4.FontSize = 13;
h4.Title = 'Br/Btot > 0.5';
h4.XLabel = '\epsilon';
h4.YLabel = '\eta';
h4.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h4.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h4.CellLabelColor = 'none';
h4.MissingDataColor = [1 1 1];
if any(isnan(h4.ColorData(:)))
    h4.MissingDataLabel = 'no resistant';
end

%% Heatmap of crossing times for Cr/Ctot (epsilon on columns, eta on rows)
nexttile(tl_row2, col2, [1 heatmap_span]);
h5 = heatmap(epsilon_labels, eta_labels, Cr_Ctot_time_1', 'Colormap', parula);
h5.FontSize = 13;
h5.Title = 'Cr/Ctot > 0.5';
h5.XLabel = '\epsilon';
h5.YLabel = '\eta';
h5.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h5.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h5.CellLabelColor = 'none';
h5.MissingDataColor = [1 1 1];
if any(isnan(h5.ColorData(:)))
    h5.MissingDataLabel = 'no resistant';
end

%% Heatmap of crossing times for Br and Cr (epsilon on columns, eta on rows)
hs1 = struct(h1);
hs1.Axes.YAxisLocation = 'left';
hs1.Axes.YLabel.Rotation = 0;
hs1.Axes.Title.FontWeight = 'normal';
ylabel(hs1.Colorbar, sprintf("Time \n when Br \n and Cr \n exceed \n 3 CFU"));
hs1.Colorbar.Label.Rotation = 0;
hs1.Colorbar.Label.HorizontalAlignment = 'center';
hs1.Colorbar.Label.Position(1) = hs1.Colorbar.Label.Position(1) + 0;
hs1.Colorbar.Label.Position(2) = hs1.Colorbar.Label.Position(2) + 15;

%% Heatmap of crossing times for Br (epsilon on columns, eta on rows)
hs2 = struct(h2);
hs2.Axes.Title.FontWeight = 'normal';
hs2.Axes.YLabel.Rotation = 0;
ylabel(hs2.Colorbar, sprintf("Time \n when Br \n exceeds \n 3 CFU"));
hs2.Colorbar.Label.Rotation = 0;
hs2.Colorbar.Label.HorizontalAlignment = 'center';
hs2.Colorbar.Label.Position(1) = hs2.Colorbar.Label.Position(1) + 0;
hs2.Colorbar.Label.Position(2) = hs2.Colorbar.Label.Position(2) + 15;

%% Heatmap of crossing times for Cr (epsilon on columns, eta on rows)
hs3 = struct(h3);
hs3.Axes.Title.FontWeight = 'normal';
hs3.Axes.YLabel.Rotation = 0;
ylabel(hs3.Colorbar, sprintf("Time \n when Cr \n exceeds \n 3 CFU"));
hs3.Colorbar.Label.Rotation = 0;
hs3.Colorbar.Label.HorizontalAlignment = 'center';
hs3.Colorbar.Label.Position(1) = hs3.Colorbar.Label.Position(1) + 0;
hs3.Colorbar.Label.Position(2) = hs3.Colorbar.Label.Position(2) + 15;

%% Heatmap of crossing times for Br/Btot (epsilon on columns, eta on rows)
hs4 = struct(h4);
hs4.Axes.YAxisLocation = 'left';
hs4.Axes.YLabel.Rotation = 0;
hs4.Axes.Title.FontWeight = 'normal';
ylabel(hs4.Colorbar, sprintf("Time \n when \n Br/Btot \n exceeds \n 0.5"));
hs4.Colorbar.Label.Rotation = 0;
hs4.Colorbar.Label.HorizontalAlignment = 'center';
hs4.Colorbar.Label.Position(1) = hs4.Colorbar.Label.Position(1) + 0;
hs4.Colorbar.Label.Position(2) = hs4.Colorbar.Label.Position(2) + 15;

%% Heatmap of crossing times for Cr/Ctot and half alpha (epsilon on columns, eta on rows)
hs5 = struct(h5);
hs5.Axes.YAxisLocation = 'left';
hs5.Axes.YLabel.Rotation = 0;
hs5.Axes.Title.FontWeight = 'normal';
ylabel(hs5.Colorbar, sprintf("Time \n when \n Cr/Ctot \n exceeds \n 0.5"));
hs5.Colorbar.Label.Rotation = 0;
hs5.Colorbar.Label.HorizontalAlignment = 'center';
hs5.Colorbar.Label.Position(1) = hs5.Colorbar.Label.Position(1)+0;
hs5.Colorbar.Label.Position(2) = hs5.Colorbar.Label.Position(2)+15;

%% Titles
h7 = title(tl_row2, 'Effect of \epsilon and \eta on half-population resistance acquisition time');
h7.FontWeight = 'normal';
h7.FontSize = 16;

h6 = sgtitle(tl, 'Effect of \epsilon and \eta on resistance acquisition time');
h6.FontSize = 18;

