%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fits Nextnano obtained Transition Energies at the simulated QW Widths
% with a power curve. QW Widths for each structure observed by PL is
% extracted by interpolating this power fitting at P = 47% for the GaAsP
% barrier. The data is then plot for creation of figure 4 in the ACS
% publication.
%%% ________________________ Order of Process _________________________ %%%
% A. Load data from nanoskived_results.h5
% B. Obtain medians and quartiles from loaded data
% C. Fit Nextnano QW Width (Y) vs Transition Energy (X) w/ Power curve
% D. Interpolate QW Widths from above fitting
% E. Prepare and Plot ACS figure 4
% (F.) Extra figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% A. Load Data
QW_vary = h5read('..\..\h5 Format\nanoskived_results.h5','/results/Nextnano/QW_Vector');
Strained_Energies = zeros(numel(QW_vary),3); % Columns: 44%, 47% and 50% P Barrier content
    Strained_Energies(:,1) = h5read('..\..\h5 Format\nanoskived_results.h5','/results/Nextnano/Strained_0.44_P_Energy');
    Strained_Energies(:,2) = h5read('..\..\h5 Format\nanoskived_results.h5','/results/Nextnano/Strained_0.47_P_Energy');
    Strained_Energies(:,3) = h5read('..\..\h5 Format\nanoskived_results.h5','/results/Nextnano/Strained_0.50_P_Energy');
Unstrained_Energies = h5read('..\..\h5 Format\nanoskived_results.h5','/results/Nextnano/Unstrained_0.47_P_Energy');

%Import master table and remove flagged data points (due to poor LSW fitting)
master_table = h5_import_pl_master_table('..\..\h5 Format\nanoskived_results.h5');
master_table = master_table(~master_table.Flag,:);

f_thick = [200,250,300,300,300,     300,300,300,300,    500,   300,300,300,300,300,300,     200]; %Heights of all flakes
f_pos = 5000 + cumsum(f_thick) - f_thick/2;
f_pos = f_pos./1e3;

%% B. PL Medians (Create table for error bars based on medians and quartiles)
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
    if j == "Height" %Skip columns that don't require error analysis
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

% Removal of Data from contaminated slices (High rejection rate or poor visibility during measurement - S1F1 and S3F8)
PL_rem = ~ismember(string(master_table.Flake),["S1F1","S3F8"]);
master_table = master_table(PL_rem,:);
PL_rem_mt = ismember(median_table.Height,f_pos([2:4,6:9,11:15,17]));
median_table = median_table(PL_rem_mt,:);

%% C. and D. QW Width Interpolation
sz = size(Strained_Energies);
EL_fits = zeros(sz(1),sz(2));

% with 44% P concentration
[xData, yData] = prepareCurveData( Strained_Energies(:,1), QW_vary.' );
% Set up fittype and options.
ft = fittype( 'power2' ); %interpolation
% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft );
L_lower = fitresult(1.546-0.003);

% with 47% P concentration
[xData, yData] = prepareCurveData( Strained_Energies(:,2), QW_vary.' );
% Set up fittype and options.
ft = fittype( 'power2' ); %interpolation
% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft );
master_table.QW_Width = fitresult(master_table.Transition_Energy);
cats = categories(master_table.Flake);

%Update QW widths by interpolating curve from nextnano [QW Width vs Transition Energy]
for i = 1:numel(median_table.Height)
    vals = master_table.QW_Width(master_table.Height == median_table.Height(i));
    medians = prctile(vals,[25 50 75]);
    median_table.QW_Width(i,:) = [medians(1)-medians(2) medians(2) medians(3)-medians(2)];
end
save("..\..\results\master_table_NN.mat","master_table","-mat");
L_mid = fitresult(1.546);

% with 50% P concentration
[xData, yData] = prepareCurveData( Strained_Energies(:,3), QW_vary.' );
% Set up fittype and options.
ft = fittype( 'power2' ); %interpolation
% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft );
L_upper = fitresult(1.546+0.004);

disp(strcat('L = ',num2str(L_mid),"(+",num2str(L_upper-L_mid),")","(-",num2str(L_mid-L_lower),")"))

%% E.
% Invert Axis ( E[Y] vs L[X] )
% with 44% P concentration
fitresult = curve_fit( Strained_Energies(:,1) , QW_vary.');
fit_lower = fitresult;

% with 47% P concentration
fitresult = curve_fit( Strained_Energies(:,2) , QW_vary.');
fit_mid = fitresult;

% with 50% P concentration
fitresult = curve_fit( Strained_Energies(:,3) , QW_vary.');
fit_upper = fitresult;
% Unstrained
fitresult = curve_fit(Unstrained_Energies, QW_vary.');
unstrained_fit = fitresult;

% Fits for Differential Fittings (Fig. 4c)
QW_Array = linspace(3,12,1000);
% with 47% P concentration
Dfit_Strained = curve_fit( QW_vary.' , Strained_Energies(:,2));
dY_Strain = differentiate(Dfit_Strained,QW_Array);
X_Strain_Energies = Dfit_Strained(QW_Array);

% unstrained
Dfit_Unstrained = curve_fit( QW_vary.' , Unstrained_Energies);
dY_Unstrain = differentiate(Dfit_Unstrained,QW_Array);
X_Unstrain_Energies = Dfit_Unstrained(QW_Array);

% Figure 4 ACS Nano _______________________________________________________
f = figure(2);
clf

% ============================== Scatter ================================ %
subplot(6,3,1:3);
colormap hsv
swarmchart(master_table,'Height','QW_Width','filled','ColorVariable','Flake',...
    'MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2,'XJitterWidth',0.2,'XJitter','randn')
hold on
x_med = median_table.Height;
vals = median_table.QW_Width;
errorbar(x_med,vals(:,2),vals(:,1),vals(:,3),'kx','MarkerSize',12,'LineWidth',2)
y = master_table.QW_Width;
x = master_table.Height;
plot_preds(x,y);
ylabel({'QW Width';'(nm)'})
ylim([7 9])
xlabel({'Slice Height (\mum)',''})
xlim([4.9 10.2])
box on
set(gca,'FontSize',18,'TickDir','in','linewidth',2)
hold off

% ============================== E vs QW ================================ %
subplot(6,3,4:12);
hold on
title({'';''})
% [RED] Patch for Energy
xline(1.546,'r-','LineWidth',2)
y = [0 12 12 0];
x = [1.543 1.543 1.55 1.55];
patch(x,y,'red','EdgeColor','none','FaceAlpha',0.2)

% [GREEN] Patch(es) for QW Width (IQR and P composition)
%1. Green
[QW_meds,QW_iqr] = quartiles(master_table.QW_Width);
yline(L_mid,'k',LineWidth=2)
y = [QW_meds(2)+QW_meds(1) QW_meds(2)+QW_meds(3) QW_meds(2)+QW_meds(3) QW_meds(2)+QW_meds(1)];
x = [0 0 2 2];
patch(x,y,'green','EdgeColor','none','FaceAlpha',0.3)
%2. Green
y = [L_lower L_upper L_upper L_lower];
x = [0 0 2 2];
patch(x,y,'green','EdgeColor','none','FaceAlpha',0.15)

plot(Strained_Energies(:,1),QW_vary,'ko','linewidth',2,'MarkerSize',8)
plot(Strained_Energies(:,2),QW_vary,'bx','linewidth',2,'MarkerSize',12)
plot(Strained_Energies(:,3),QW_vary,'kd','linewidth',2,'MarkerSize',8)
plot(fit_lower, 'k--')
plot(fit_mid, 'k')
plot(fit_upper, 'k--')

hold off
a=get(gca,'Children'); %obtain data tags
legend([a(11),a(6),a(5),a(4)],{'Experiment (1.546 eV)','44% P','47% P','50% P'}, ...
    'location','best','FontSize',18)
xlabel('Energy (eV)')
ylabel({'QW Thickness','(nm)'})
axis tight
xlim([1.48 1.62])
ylim([6 11])
box on
set(gca,'FontSize',18,'TickDir','in','linewidth',2)

% =========================== Differential ============================== %
subplot(6,3,13:18);
hold on
plot(1e3*X_Strain_Energies-1546,1e3*dY_Strain,LineWidth=2)
plot(1e3*X_Unstrain_Energies-1546,1e3*dY_Unstrain,LineWidth=2)
% [RED] Patch for Energy
xline(0,'r-','LineWidth',2)
y = [-100 0 0 -100];
x = [-3 -3 4 4];
patch(x,y,'red','EdgeColor','none','FaceAlpha',0.2)
hold off
b=get(gca,'Children'); %obtain data tags
 legend([b(4),b(3)],{'Strained','Unstrained'}, ...
     'location','best','FontSize',18)
xlabel('Energy wrt. 1546 meV (meV)')
ylabel({'dE/dL','(meV/nm)'})
xlim([-30 30])
ylim([-60 -17])
box on
set(gca,'FontSize',18,'TickDir','in','linewidth',2)
movegui southeast

%% (F.) Strained vs Unstrained
figure(3)
clf
hold on
plot(QW_vary.',Strained_Energies(:,2),'bo','linewidth',2)
plot(QW_vary.',Unstrained_Energies,'ro','linewidth',2)
plot(Dfit_Strained)%,'b-.','LineWidth',2)
plot(Dfit_Unstrained)%,'r-.','LineWidth',2)
hold off
c=get(gca,'Children'); %obtain data tags
legend([c(4),c(3)],{'Strained','Unstrained'}, ...
    'location','best','FontSize',18)
ylabel('Energy (eV)')
xlabel('QW Thickness (nm)')
axis tight
box on
set(gca,'FontSize',18,'TickDir','in','linewidth',2)

%% Functions
function g_bounds = plot_preds(x,y)
    fitresult = fit(x,y,'poly1');
    ci = confint(fitresult,0.95);
    plot(x,(fitresult.p1).*x+fitresult.p2,'k','LineWidth',2);

    g_bounds = [ci(1,1)-fitresult.p1,fitresult.p1,ci(2,1)-fitresult.p1];
end

function fitresult = curve_fit(xData,yData) %poly fit
    [xData, yData] = prepareCurveData( xData, yData );
    ft = fittype( 'power2' );
    [fitresult, ~] = fit( xData, yData, ft );
end