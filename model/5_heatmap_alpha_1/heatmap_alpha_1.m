%% Heatmap realization for third drug community
% In this file a heatmap is created with different values of epsilon and
% eta. Alpha (the amino acid production rate) is increased and decreased to
% show changes over time in resistance acquisition. The file may take 
% several minutes to process. To reduce the running time, alpha = alpha * 1.5 
% can be reduced to alpha = alpha * 1.2.

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
Br_Btot_time_2  = nan(n,n);
Br_Btot_time_3  = nan(n,n);
Cr_Ctot_time_1  = nan(n,n);
Cr_Ctot_time_2  = nan(n,n);
Cr_Ctot_time_3  = nan(n,n);

for i = 1:n
    epsilon = epsilon_values(i);
    for k = 1:n
        eta = eta_values(k);

        time = cfu(:,1);
        cfu_A = cfu(:, 11);
        cfu_B = cfu(:, 12);
        cfu_C = cfu(:, 13);
        
        for j = 1:3
            %% Parameters from no-drug communities of ampicillin experiments
            mu_A = 0.0003;
            mu_Ar = mu_A * scale_factor;
            mu_B = 0.0005;
            mu_Br = mu_B * scale_factor;
            mu_C = 0.0004;
            mu_Cr = mu_C * scale_factor;
            switch j
                case 1
                    alpha_A = 0.7645/2;
                    alpha_B = 0.0567/2;
                    alpha_C = 0.1202/2;
                case 2
                    alpha_A = 0.7645;
                    alpha_B = 0.0567;
                    alpha_C = 0.1202;
                case 3
                    alpha_A = 0.7645*1.5;
                    alpha_B = 0.0567*1.5;
                    alpha_C = 0.1202*1.5;
            end
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
            switch j
                case 1
                    idx_cross = find(Br > 1 & Cr > 1, 1, 'first');
                case {2,3}
                    idx_cross = find(Br > 3 & Cr > 3, 1, 'first');
            end
            idx_Br    = find(Br_Btot > 0.5, 1, 'first');
            idx_Cr    = find(Cr_Ctot > 0.5, 1, 'first');
            
            switch j
                case 1
                    if ~isempty(idx_cross), crossing_time_1(i,k) = t(idx_cross); end
                    if ~isempty(idx_Br),    Br_Btot_time_1(i,k)  = t(idx_Br);    end
                    if ~isempty(idx_Cr),    Cr_Ctot_time_1(i,k)  = t(idx_Cr);    end
                case 2
                    if ~isempty(idx_cross), crossing_time_2(i,k) = t(idx_cross); end
                    if ~isempty(idx_Br),    Br_Btot_time_2(i,k)  = t(idx_Br);    end
                    if ~isempty(idx_Cr),    Cr_Ctot_time_2(i,k)  = t(idx_Cr);    end
                case 3
                    if ~isempty(idx_cross), crossing_time_3(i,k) = t(idx_cross); end
                    if ~isempty(idx_Br),    Br_Btot_time_3(i,k)  = t(idx_Br);    end
                    if ~isempty(idx_Cr),    Cr_Ctot_time_3(i,k)  = t(idx_Cr);    end
            end
        end
    end

end

%% Heatmap plot (epsilon on columns, eta on rows)
epsilon_labels = arrayfun(@(x) sprintf('%.4g', x), epsilon_values, 'UniformOutput', false);
eta_labels = arrayfun(@(x) sprintf('%.4g', x), eta_values, 'UniformOutput', false);

figure;
tl = tiledlayout(3,3, 'TileSpacing','tight','Padding','compact');
tl.OuterPosition = [0 0 0.94 1];

%% Heatmap of crossing times for half alpha (epsilon on columns, eta on rows)
nexttile;
h1 = heatmap(epsilon_labels, eta_labels, crossing_time_1', 'Colormap', parula);
h1.FontSize = 13;
h1.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h1.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h1.Title = '\alpha/2';
h1.XLabel = '\epsilon';
h1.YLabel = '\eta';
h1.CellLabelColor = 'none';
h1.MissingDataColor = [1 1 1];
if any(isnan(h1.ColorData(:)))     
    h1.MissingDataLabel = 'no resistant'; 
end

%% Heatmap of crossing times for normal alpha (epsilon on columns, eta on rows)
nexttile;
h2 = heatmap(epsilon_labels, eta_labels, crossing_time_2', 'Colormap', parula);
h2.FontSize = 13;
h2.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h2.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h2.Title = '\alpha';
h2.XLabel = '\epsilon';
h2.YLabel = '\eta';
h2.CellLabelColor = 'none';
h2.MissingDataColor = [1 1 1];
if any(isnan(h2.ColorData(:)))     
    h2.MissingDataLabel = 'no resistant'; 
end

%% Heatmap of crossing times for double alpha (epsilon on columns, eta on rows)
nexttile;
h3 = heatmap(epsilon_labels, eta_labels, crossing_time_3', 'Colormap', parula);
h3.FontSize = 13;
h3.Title = '\alpha \times 1.5';
h3.XLabel = '\epsilon';
h3.YLabel = '\eta';
h3.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h3.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h3.CellLabelColor = 'none';
h3.MissingDataColor = [1 1 1];
if any(isnan(h3.ColorData(:)))     
    h3.MissingDataLabel = 'no resistant'; 
end

%% Heatmap of crossing times for Br/Btot and half alpha (epsilon on columns, eta on rows)
nexttile;
h4 = heatmap(epsilon_labels, eta_labels, Br_Btot_time_1', 'Colormap', parula);
h4.FontSize = 13;
h4.XLabel = '\epsilon';
h4.YLabel = '\eta';
h4.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h4.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h4.CellLabelColor = 'none';
h4.MissingDataColor = [1 1 1];
if any(isnan(h4.ColorData(:)))     
    h4.MissingDataLabel = 'no resistant'; 
end

%% Heatmap of crossing times for Br/Btot and normal alpha (epsilon on columns, eta on rows)
nexttile;
h5 = heatmap(epsilon_labels, eta_labels, Br_Btot_time_2', 'Colormap', parula);
h5.FontSize = 13;
h5.XLabel = '\epsilon';
h5.YLabel = '\eta';
h5.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h5.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h5.CellLabelColor = 'none';
h5.MissingDataColor = [1 1 1];
if any(isnan(h5.ColorData(:)))     
    h5.MissingDataLabel = 'no resistant'; 
end

%% Heatmap of crossing times for Br/Btot and double alpha (epsilon on columns, eta on rows)
nexttile;
h6 = heatmap(epsilon_labels, eta_labels, Br_Btot_time_3', 'Colormap', parula);
h6.FontSize = 13;
h6.XLabel = '\epsilon';
h6.YLabel = '\eta';
h6.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h6.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h6.CellLabelColor = 'none';
h6.MissingDataColor = [1 1 1];
if any(isnan(h6.ColorData(:)))     
    h6.MissingDataLabel = 'no resistant'; 
end

%% Heatmap of crossing times for Cr/Ctot and half alpha (epsilon on columns, eta on rows)
nexttile;
h7 = heatmap(epsilon_labels, eta_labels, Cr_Ctot_time_1', 'Colormap', parula);
h7.FontSize = 13;
h7.XLabel = '\epsilon';
h7.YLabel = '\eta';
h7.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h7.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h7.CellLabelColor = 'none';
h7.MissingDataColor = [1 1 1];
if any(isnan(h7.ColorData(:)))     
    h7.MissingDataLabel = 'no resistant'; 
end

%% Heatmap of crossing times for Cr/Ctot and normal alpha(epsilon on columns, eta on rows)
nexttile;
h8 = heatmap(epsilon_labels, eta_labels, Cr_Ctot_time_2', 'Colormap', parula);
h8.FontSize = 13;
h8.XLabel = '\epsilon';
h8.YLabel = '\eta';
h8.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h8.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h8.CellLabelColor = 'none';
h8.MissingDataColor = [1 1 1];
if any(isnan(h8.ColorData(:)))     
    h8.MissingDataLabel = 'no resistant'; 
end

%% Heatmap of crossing times for Cr/Ctot and double alpha(epsilon on columns, eta on rows)
nexttile;
h9 = heatmap(epsilon_labels, eta_labels, Cr_Ctot_time_3', 'Colormap', parula);
h9.FontSize = 13;
h9.XLabel = '\epsilon';
h9.YLabel = '\eta';
h9.XDisplayLabels = repmat({''}, 1, numel(epsilon_labels));
h9.YDisplayLabels = repmat({''}, 1, numel(eta_values));
h9.CellLabelColor = 'none';
h9.MissingDataColor = [1 1 1];
if any(isnan(h9.ColorData(:)))     
    h9.MissingDataLabel = 'no resistant'; 
end

%% Heatmap of crossing times for half alpha (epsilon on columns, eta on rows)
hs1 = struct(h1);
hs1.Axes.YAxisLocation = 'left';
hs1.Axes.YLabel.String = 'Br & Cr > 3 CFU';
hs1.Axes.YLabel.Rotation = 0;
hs1.Axes.Title.FontWeight = 'normal';

%% Heatmap of crossing times for normal alpha (epsilon on columns, eta on rows)
hs2 = struct(h2);
hs2.Axes.Title.FontWeight = 'normal';
hs2.Axes.YLabel.Rotation = 0;

%% Heatmap of crossing times for double alpha (epsilon on columns, eta on rows)
hs3 = struct(h3);
hs3.Axes.Title.FontWeight = 'normal';
hs3.Axes.YLabel.Rotation = 0;
ylabel(hs3.Colorbar, sprintf("Time \n when Br \n and Cr \n exceed \n 3 CFU"));
hs3.Colorbar.Label.Rotation = 0;
hs3.Colorbar.Label.HorizontalAlignment = 'center';
hs3.Colorbar.Label.Position(1) = hs3.Colorbar.Label.Position(1) + 0;
hs3.Colorbar.Label.Position(2) = hs3.Colorbar.Label.Position(2) + 10;

%% Heatmap of crossing times for Br/Btot and half alpha (epsilon on columns, eta on rows)
hs4 = struct(h4);
hs4.Axes.YAxisLocation = 'left';
hs4.Axes.YLabel.String = 'Br/Btot > 0.5';
hs4.Axes.YLabel.Rotation = 0;

%% Heatmap of crossing times for Br/Btot and normal alpha (epsilon on columns, eta on rows)
hs5 = struct(h5);
hs5.Axes.YLabel.Rotation = 0;

%% Heatmap of crossing times for Br/Btot and double alpha (epsilon on columns, eta on rows)
hs6 = struct(h6);
hs6.Axes.YLabel.Rotation = 0;
ylabel(hs6.Colorbar, sprintf("Time \n when \n Br/Btot \n exceeds \n 0.5"));
hs6.Colorbar.Label.Rotation = 0;
hs6.Colorbar.Label.HorizontalAlignment = 'center';
hs6.Colorbar.Label.Position(1) = hs6.Colorbar.Label.Position(1) + 0;
hs6.Colorbar.Label.Position(2) = hs6.Colorbar.Label.Position(2) + 15;

%% Heatmap of crossing times for Cr/Ctot and half alpha (epsilon on columns, eta on rows)
hs7 = struct(h7);
hs7.Axes.YAxisLocation = 'left';
hs7.Axes.YLabel.String = 'Cr/Ctot > 0.5';  
hs7.Axes.YLabel.Rotation = 0;

%% Heatmap of crossing times for Cr/Ctot and normal alpha(epsilon on columns, eta on rows)
hs8 = struct(h8);
hs8.Axes.YLabel.Rotation = 0;

%% Heatmap of crossing times for Cr/Ctot and double alpha (epsilon on columns, eta on rows)
hs9 = struct(h9);
hs9.Axes.YLabel.Rotation = 0;
ylabel(hs9.Colorbar, sprintf("Time \n when \n Cr/Ctot \n exceeds \n 0.5"));
hs9.Colorbar.Label.Rotation = 0;
hs9.Colorbar.Label.HorizontalAlignment = 'center';
hs9.Colorbar.Label.Position(1) = hs9.Colorbar.Label.Position(1)+0;
hs9.Colorbar.Label.Position(2) = hs9.Colorbar.Label.Position(2)+15;

h10 = sgtitle('Effect of amino acid production rate on resistance acquisition time');
h10.FontSize = 18;

