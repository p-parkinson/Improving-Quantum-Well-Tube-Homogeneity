%% Initialise Tables
n=height(cell2mat(QWs_SuperTable.Region_Size(:,1)));
master_celltable = table('Size', [n,17],'VariableTypes',...
    ["string","double", "double","double",...
    "table","double","double","double",...
    "logical","double","double","double","double","double" ...
    "double","double","double"]);

master_celltable.Properties.VariableNames = ["Flake","Height","FlakeID","Size",...
    "XY","PL_Data","LSW_Fit","RMS",...
    "Flag","QW_Width","Transition_Energy","Urbach_Energy","Peak_Emission","Fermi_Energy",...
    "Sigma", "Temperature", "Intensity"];

%Store XY Co-ordinates of Individual Regions
xy_table = table('Size',[n,2],'VariableTypes',["double","double"]);
xy_table.Properties.VariableNames = ["X","Y"];

PLs_table = cell(n,1); %PL spectra table
PLf_table = cell(n,1); %PL LSW fitting table
PLr_table = cell(n,1); %PL RMS table

%% Populate Tables
k=1;
for i =1:height(QWs_SuperTable)
    mask = QWs_SuperTable{i,'Flagged'}{1}; %if flagged with (1), then this is a rejected fitting
    if i == 5
        continue
    end

    no_structures = numel(mask);

    % Metadata
    master_celltable{k:k+no_structures-1,"Flake"} = QWs_SuperTable{i,'Flake'};
    master_celltable{k:k+no_structures-1,"Height"} = QWs_SuperTable{i,'Height'}{1};
    master_celltable{k:k+no_structures-1,"FlakeID"} = QWs_SuperTable{i,'Region_ID'}{1};

    % Interwire and Intrawire Data
    master_celltable{k:k+no_structures-1,"Size"} = QWs_SuperTable{i,'Region_Size'}{1};
    for j=1:2
        xy_table{k:k+no_structures-1,j} = QWs_SuperTable{i,'XYZ'}{1}(:,j);
    end
    
    % PL Spectra and Fitting Data
    ctr = 1;
    for w = k:k+no_structures-1
        temp = QWs_SuperTable{i,'PL_Spectra'}{:}(:,:);
        PLs_table(w,:) = {temp(ctr,:)};

        temp = QWs_SuperTable{i,'PL_Fit'}{:}(:,:);
        PLf_table(w,:) = {temp(ctr,:)};

        temp = QWs_SuperTable{i,'RMS'}{:}(:,:);
        PLr_table(w,:) = {temp(ctr,:)};

        ctr = ctr + 1;
    end
    
    % LSW Fitting Outputs
    master_celltable{k:k+no_structures-1,"Flag"} = logical(QWs_SuperTable{i,'Flagged'}{1});
    master_celltable{k:k+no_structures-1,"QW_Width"} = QWs_SuperTable{i,'QW_Width'}{1}(:,1);  
    master_celltable{k:k+no_structures-1,"Transition_Energy"} = QWs_SuperTable{i,'Transition_Energy'}{1}(:,1);
    master_celltable{k:k+no_structures-1,"Urbach_Energy"} = QWs_SuperTable{i,'Urbach_Energy'}{1}(:,1).*1e3;
    master_celltable{k:k+no_structures-1,"Peak_Emission"} = QWs_SuperTable{i,'Wavelength'}{1}(:,1);
    master_celltable{k:k+no_structures-1,"Fermi_Energy"} = QWs_SuperTable{i,'Fermi_Energy'}{1}(:,1);
    master_celltable{k:k+no_structures-1,"Intensity"} = QWs_SuperTable{i,'Intensity'}{1}(:,1);
    master_celltable{k:k+no_structures-1,"Sigma"} = QWs_SuperTable{i,'Sigma'}{1}(:,1).*1e3;
    master_celltable{k:k+no_structures-1,"Temperature"} = QWs_SuperTable{i,'Temperature'}{1}(:,1);

    k = k+no_structures;
end

%% Post-Process Populated Tables (To correct datatype)
master_celltable.Flake = categorical(master_celltable.Flake);
master_celltable.XY = xy_table;

master_celltable.PL_Data = PLs_table;
master_celltable.LSW_Fit = PLf_table;
master_celltable.RMS = PLr_table;

master_table = table(master_celltable{:,"Flake"},...
    master_celltable{:,"Height"},...
    master_celltable{:,"FlakeID"},...
    master_celltable{:,"Size"},...
    table2array(master_celltable{:,"XY"}),...
    master_celltable{:,"PL_Data"},...
    master_celltable{:,"LSW_Fit"},...
    cell2mat(master_celltable{:,"RMS"}),...
    master_celltable{:,"Flag"},...
    master_celltable{:,"QW_Width"},...
    master_celltable{:,"Transition_Energy"},...
    master_celltable{:,"Urbach_Energy"},...
    master_celltable{:,"Peak_Emission"},...
    master_celltable{:,"Fermi_Energy"},...
    master_celltable{:,"Sigma"},...
    master_celltable{:,"Temperature"},...
    master_celltable{:,"Intensity"},...
    'VariableNames',["Flake","Height","FlakeID","Size",...
    "XY",...
    "PL_Data","LSW_Fit","RMS",...
    "Flag","QW_Width","Transition_Energy","Urbach_Energy","Peak_Emission","Fermi_Energy",...
    "Sigma", "Temperature", "Intensity"]);

% Save Table
save("..\..\results\master_table_LSW.mat","master_table");