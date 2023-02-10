%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Swarm plots of parameters obtained from LSW fitting of PL data along
% with morphological data from SEM image analysis. Includes plots used in
% figures 2. and 3. in the ACS publication.
%%% ________________________ Order of Process _________________________ %%%
% A. Load data from nanoskived_results.h5
% B. Obtain medians and quartiles from loaded data
% C. Swarm plots (Figure 3 ACS)
% D. Transition energy histogram (Figure 2 ACS)
% E. Example PL curve (Figure 2 ACS)
% (F.) Extra swarm plots
% G. Errors on slopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialisation
f_thick = [200,250,300,300,300,     300,300,300,300,    500,   300,300,300,300,300,300,     200]; %Heights of all flakes
f_pos = 5000 + cumsum(f_thick) - f_thick/2;
f_pos = f_pos./1e3;

lns = hsv(17); %colormap for swarm plots

%% A. Load Data Tables
% PL Data Master Table
master_table = h5_import_pl_master_table('..\..\h5 Format\nanoskived_results.h5');
master_table = master_table(~master_table.Flag,:);

% SEM Data Master Table
SEM_master = h5_import_sem_master_table('..\..\h5 Format\nanoskived_results.h5');

%% B. Calculate Medians
% PL Medians
% Obtain medians for each data variable
median_table = table('Size', [numel(f_pos([1:4,6:9,11:17])),7],'VariableTypes',...
    ["double", ...
    "table","table", ...
    "table","table", ...
    "table","table"]);
median_table.Properties.VariableNames = ["Height", ...
    "QW_Width","Transition_Energy",...
    "Urbach_Energy","Sigma",...
    "Temperature","Intensity"];

median_table.("Height") = f_pos([1:4,6:9,11:17]).'; %No PL in (S1F1), S1F5, TEM
for j = string(median_table.Properties.VariableNames)
    T = zeros(height(median_table),3);
    if j == "Height"
        continue
    end
    k = 1;
    for i = [1:4,6:9,11:17]
        Tset = master_table(master_table.Height == f_pos(i),:);
        [Qs,~] = quartiles(Tset.(j)(~Tset.Flag,:));
        T(k,:) = Qs;
        k = k + 1;
    end
    median_table.(j) = T;
end

% SEM Medians
% Obtain medians for each data variable
flakes = string(categories(SEM_master.Flake));
flakes = flakes(2:end-1); %Remove "N/A" at (1) and "sample" at (end)

SEM_medians = table('Size', [numel(f_pos([1:4,6:9,11:17])),3],'VariableTypes',...
    ["double","double","double"]);
SEM_medians.Properties.VariableNames = ["Height","Neighbour Distance",...
    "Diameters"];
SEM_medians.Height = f_pos([1:4,6:9,11:17]).'; %No SEM for (S1F1), S1F5, and TEM
NNs = zeros(numel(SEM_medians.Height),3);
diams = zeros(numel(SEM_medians.Height),3);
diam_IQR = zeros(numel(SEM_medians.Height),1);

k = 1; %counter
for i = [1:4,6:9,11:17]
    T = SEM_master.("Nearest Neighbour");
    T = T(SEM_master.Height == f_pos(i));
    QRTs = quantile(T,[0.25 0.50 0.75]);
    NNs(k,:) = [QRTs(1)-QRTs(2) QRTs(2) QRTs(3)-QRTs(2)];

    T = SEM_master.("Diameters");
    T = T(SEM_master.Height == f_pos(i));
    QRTs = quantile(T,[0.25 0.50 0.75]);
    diams(k,:) = [QRTs(1)-QRTs(2) QRTs(2) QRTs(3)-QRTs(2)];
    diam_IQR(k) = iqr(T);

    k = k + 1;
end
SEM_medians.("Neighbour Distance") = NNs;
SEM_medians.("Diameters") = diams;

% Removal of Data from contaminated slices (High rejection rate or poor visibility during measurement)
PL_rem = ~ismember(string(master_table.Flake),["S1F1","S3F8"]);
master_table = master_table(PL_rem,:);
PL_rem_mt = ismember(median_table.Height,f_pos([2:4,6:9,11:15,17]));
median_table = median_table(PL_rem_mt,:);

SEM_rem = ~ismember(string(SEM_master.Flake),"S1F5");
SEM_master = SEM_master(SEM_rem,:);
SEM_rem_mt = ismember(SEM_medians.Height,f_pos([1:4,6:9,11:17]));
SEM_medians = SEM_medians(SEM_rem_mt,:);
    %% C. Swarm Plots
    
    %%% PL Violins
    f = figure(1);
    clf
    tiledlayout(6,1)
    median_iqr = zeros(6,2); %[median_iqr}
    bounds = zeros(7,3);

    % EMISSION
    ax = nexttile; colormap(lns); y_lims = [1.53 1.57]; num_ys = 3;
    swarmchart(master_table,'Height','Transition_Energy','filled','ColorVariable','Flake',...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'XJitterWidth',0.2,'XJitter','randn')
    hold on
    vals = median_table.Transition_Energy;
    x_med = median_table.Height;
    errorbar(x_med,vals(:,2),vals(:,1),vals(:,3),'kx','MarkerSize',9,'LineWidth',2)
    y = master_table.Transition_Energy;
    x = master_table.Height;
    median_iqr(1,:) = [median(y),iqr(y)];
    bounds(1,:) = plot_preds(x,y);
    hold off
    ylabel({'Transition';'Energy';'(eV)'})
    ylim(y_lims); yticks(ax,linspace(y_lims(1),y_lims(2),num_ys))
    xlabel('');xlim([4.9 10.2])
    box on
    yyaxis right
    ylabel('PL')
    set(gca,'FontSize',18,'TickDir','in','linewidth',2,'XTickLabel','','YTick',[1.52,1.54,1.56],'YTickLabel',[1.52,1.54,1.56],'YColor','k')
    
    % TEMPERATURE
    ax = nexttile; colormap(lns); y_lims = [280 450]; num_ys = 3;
    swarmchart(master_table,'Height','Temperature','filled','ColorVariable','Flake',...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'XJitterWidth',0.2,'XJitter','randn')
    hold on
    vals = median_table.Temperature;
    x_med = median_table.Height;
    errorbar(x_med,vals(:,2),vals(:,1),vals(:,3),'kx','MarkerSize',9,'LineWidth',2)
    y = master_table.Temperature;
    x = master_table.Height;
    median_iqr(2,:) = [median(y),iqr(y)];
    bounds(2,:) = plot_preds(x,y);
    hold off
    ylabel({'Emission';'Temperature';'(K)'})
    ylim(y_lims); yticks(ax,linspace(y_lims(1),y_lims(2),num_ys))
    xlabel('');xlim([4.9 10.2])
    box on
    yyaxis right
    ylabel('PL')
    set(gca,'FontSize',18,'TickDir','in','linewidth',2,'XTickLabel','','YTick',[], 'YTickLabel', [],'YColor','k')
    
    % URBACH ENERGY
    ax = nexttile; colormap(lns); y_lims = [0 50]; num_ys = 3;
    swarmchart(master_table,'Height','Urbach_Energy','filled','ColorVariable','Flake',...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'XJitterWidth',0.2,'XJitter','randn')
    hold on
    vals = median_table.Urbach_Energy;
    x_med = median_table.Height;
    errorbar(x_med,vals(:,2),vals(:,1),vals(:,3),'kx','MarkerSize',9,'LineWidth',2)
    y = master_table.Urbach_Energy;
    x = master_table.Height;
    median_iqr(3,:) = [median(y),iqr(y)];
    bounds(3,:) = plot_preds(x,y);
    hold off
    xlabel('');xlim([4.9 10.2])
    ylim(y_lims); yticks(ax,linspace(y_lims(1),y_lims(2),num_ys))
    ylabel({'Urbach';'Energy';'(meV)'})
    box on
    yyaxis right
    ylabel('PL')
    set(gca,'FontSize',18,'TickDir','in','linewidth',2,'XTickLabel','','YTick',[], 'YTickLabel', [],'YColor','k')
    
    % Sigma
    ax = nexttile; colormap(lns); y_lims = [0 30]; num_ys = 3;
    swarmchart(master_table,'Height','Sigma','filled','ColorVariable','Flake',...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'XJitterWidth',0.2,'XJitter','randn')
    hold on
    vals = median_table.Sigma;
    x_med = median_table.Height;
    errorbar(x_med,vals(:,2),vals(:,1),vals(:,3),'kx','MarkerSize',9,'LineWidth',2)
    y = master_table.Sigma;
    x = master_table.Height;
    median_iqr(4,:) = [median(y),iqr(y)];
    bounds(4,:) = plot_preds(x,y);
    hold off
    ylabel({'Sigma';'(meV)'})
    ylim(y_lims); yticks(ax,linspace(y_lims(1),y_lims(2),4))
    xlabel('');xlim([4.9 10.2])
    box on
    yyaxis right
    ylim([0 0.02])
    ylabel('PL')
    set(gca,'FontSize',18,'TickDir','in','linewidth',2,'XTickLabel','','YTick',[], 'YTickLabel', [],'YColor','k')

    %%% SEM Violins
    % DIAMETERS
    ax = nexttile; colormap(lns); y_lims = [200 900]; num_ys = 3;
    swarmchart(SEM_master,'Height','Diameters','filled','ColorVariable','Flake',...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'XJitterWidth',0.2,'XJitter','randn')
    hold on
    y = SEM_master.("Diameters");
    x = SEM_master.Height;
    vals = SEM_medians.Diameters;
    x_med = SEM_medians.Height;
    errorbar(x_med,vals(:,2),vals(:,1),vals(:,3),'kx','MarkerSize',9,'LineWidth',2)
    bounds(5,:) = plot_preds(x,y);
    ylabel({'Diameter';'(nm)'})
    xlabel('');xlim([4.9 10.2])
    ylim(y_lims); yticks(ax,linspace(y_lims(1),y_lims(2),num_ys))
    box on
    yyaxis right
    ylabel('SEM')
    set(gca,'FontSize',18,'TickDir','in','linewidth',2,'XTickLabel','','YTick',[], 'YTickLabel', [],'YColor','k')
    
    % SEPARATION
    ax = nexttile; colormap(lns); y_lims = [0 2.5]; num_ys = 3;
    swarmchart(SEM_master,'Height','Nearest Neighbour','filled','ColorVariable','Flake',...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'XJitterWidth',0.2,'XJitter','randn')
    hold on
    y = SEM_master.("Nearest Neighbour");
    x = SEM_master.Height;
    vals = SEM_medians.("Neighbour Distance");
    x_med = SEM_medians.Height;
    errorbar(x_med,vals(:,2),vals(:,1),vals(:,3),'kx','MarkerSize',9,'LineWidth',2)
    bounds(6,:) = plot_preds(x,y);
    ylim(y_lims); yticks(ax,linspace(y_lims(1),y_lims(2),num_ys))
    ylabel({'Interwire';'Separation';'(\mu{m})'})
    xlabel('Slice Height (\mum)');xlim([4.9 10.2])
    box on
    yyaxis right
    ylabel('SEM')
    set(gca,'FontSize',18,'TickDir','in','linewidth',2,'YTick',[], 'YTickLabel', [],'YColor','k') 
    movegui southeast
    f.Position = [0 0 2*380 2*680];

%% D. EV HISTFIT FOR ABSTRACT
figure(2)
sedges = 1.53:0.002:1.57;
h = histfit(master_table.Transition_Energy,numel(sedges),'kernel');

pd = fitdist(master_table.Transition_Energy,'kernel');
h(1).FaceColor = [0.9290 0.6940 0.1250];

histogram(master_table.Transition_Energy,sedges,'Normalization','probability')
ylabel('Probability Density')
xlabel('Photon Energy (eV)')
axis tight
set(gca,'FontSize',18,'TickDir','in','linewidth',3)

%% E. PL Spectra Fitting
entry = 456; %pick a spectra to plot from master_table.mat
f = figure(4);

x = h5read('..\..\h5 Format\nanoskived_results.h5','/results/PL_table/pl_spectrum_x_axis');

% Constants
h = 6.626*10^(-34);
c = 2.998*10^(8);
eV = 1.602*10^(-19);
k_B = 8.617*10^(-5); %eV/K

% Convert wavelength (nm) x axis to energy domain (eV)
E = h*c./(eV*x.*10^(-9));

plot(E,cell2mat(master_table.PL_Data(entry)),'kx','MarkerSize',8)
hold on
plot(E,cell2mat(master_table.LSW_Fit(entry)),'Linewidth',3)
hold off
ylabel('PL Intensity (Photons, a.u.)')
xlabel('Photon Energy (eV)')
axis tight
legend(["PL Data","LSW Fit"])
set(gca,'FontSize',18,'TickDir','in','linewidth',3)


%% F.Extra Plots
    % Transition Energy w/ NEXTNANO
    figure(5)
    clf
    colormap(lns);
    swarmchart(master_table,'Height','Transition_Energy','filled','ColorVariable','Flake',...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'XJitterWidth',0.2,'XJitter','randn')
    hold on
    x_med = median_table.Height;
    vals = median_table.Transition_Energy;
    errorbar(x_med,vals(:,2),vals(:,1),vals(:,3),'kx','MarkerSize',18,'LineWidth',3)
    y = master_table.Transition_Energy;
    x = master_table.Height;
    median_iqr(1,:) = [median(y),iqr(y)];
    bounds(1,:) = plot_preds(x,y);
    hold off
    ylabel('Transition Energy (eV)')
    ylim([1.53 1.57])
    xlabel('Slice Height (\mum)');xlim([4.9 10.2])
    box on
        set(gca,'FontName','Helvetica','FontSize',24,'FontWeight','bold','TickDir','in','linewidth',3,'YColor','k')

    % Intensity plot
    figure(6)
    clf
    colormap(lns);
    swarmchart(master_table,'Height','Intensity','filled','ColorVariable','Flake',...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'XJitterWidth',0.2,'XJitter','randn')
    hold on
    x_med = median_table.Height;
    vals = median_table.Intensity;
    errorbar(x_med,vals(:,2),vals(:,1),vals(:,3),'kx','MarkerSize',9,'LineWidth',2)
    y = master_table.Intensity;
    x = master_table.Height;
    median_iqr(1,:) = [median(y),iqr(y)];
    plot_preds(x,y)
    hold off
    ylabel({'Intensity';'(Photons, arb. u)'})
    ylim([0.015 0.03])
    xlabel('Slice Height (\mum)')
    xlim([4.9 10.2])
    box on
    yyaxis right
    ylabel('PL')
    set(gca,'FontSize',18,'TickDir','in','linewidth',2,'XTickLabel','','YTick',[], 'YTickLabel', [],'YColor','k')

    % QW Width w/ NEXTNANO
    figure(7)
    clf
    colormap(lns);
    swarmchart(master_table,'Height','QW_Width','filled','ColorVariable','Flake',...
        'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'XJitterWidth',0.2,'XJitter','randn')
    hold on
    x_med = median_table.Height;
    vals = median_table.QW_Width;
    errorbar(x_med,vals(:,2),vals(:,1),vals(:,3),'kx','MarkerSize',18,'LineWidth',3)
    y = master_table.QW_Width;
    x = master_table.Height;
    median_iqr(1,:) = [median(y),iqr(y)];
    bounds(7,:) = plot_preds(x,y);
    hold off
    ylabel('QW Width (nm)')
    ylim([7 9])
    xlabel('Slice Height (\mum)');xlim([4.9 10.2])
    box on
        set(gca,'FontName','Helvetica','FontSize',24,'FontWeight','bold','TickDir','in','linewidth',3,'YColor','k')
%% G. Errors on slopes for each parameter plot as a function of NW height
slope_table = table(bounds(1,:), ...
    bounds(2,:), ...
    bounds(3,:), ...
    bounds(4,:), ...
    bounds(7,:), ...
    bounds(5,:), ...
    bounds(6,:), ...
    'VariableNames',["Emission","Temperature","Urbach Energy","Sigma","QW Width","Diameter","Interwire Separation"]);

params = ["Transition_Energy","Temperature","Urbach_Energy","Sigma","QW_Width","Diameters","Nearest Neighbour"];
med_bounds = zeros(7,3);
for i = 1:7
    if i < 6
        prc = prctile(master_table.(params(i)),[25 50 75]); med_bounds(i,:) = [prc(1)-prc(2) prc(2) prc(3)-prc(2)];
    else
        prc = prctile(SEM_master.(params(i)),[25 50 75]); med_bounds(i,:) = [prc(1)-prc(2) prc(2) prc(3)-prc(2)];
    end
end

% Medians across whole population
TOTALmedian_table = table(med_bounds(1,:), ...
    med_bounds(2,:), ...
    med_bounds(3,:), ...
    med_bounds(4,:), ...
    med_bounds(5,:), ...
    med_bounds(6,:), ...
    med_bounds(7,:), ...
    'VariableNames',["Emission (eV)","Temperature (K)", ...
    "Urbach Energy (meV)","Sigma (meV)","QW Width (nm)","Diameter (nm)","Interwire Separation (\mum)"]);

%% FUNCTIONS
function g_bounds = plot_preds(x,y) %fit slope and output confidence bounds [lower slope upper]
    fitresult = fit(x,y,'poly1');
    ci = confint(fitresult,0.95);
    plot(x,(fitresult.p1).*x+fitresult.p2,'k','LineWidth',2);

    g_bounds = [ci(1,1)-fitresult.p1, fitresult.p1, ci(2,1)-fitresult.p1];
end