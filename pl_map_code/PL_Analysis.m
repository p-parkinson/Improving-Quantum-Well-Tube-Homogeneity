%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves processed data as structures and tables to a new file for quick 
% access for plotting and analysis
% Terms "Slice"/"Flake" are used interchangeably.
%%% ________________________ Order of Process _________________________ %%%
% A. Determine Watershed                            (PL_watershed.m)
% B. Fit PL Spectra                                 (PL_LSWfit.m)
% C. Flag Unphysical Fits and Associated Regions    (Section 4.)
% D. Calculate Co-ordinates of Regions              (PL_coordinates.m)
% E. Save Data in Structures and Tables
%
% ---> Then send to "table_processing.m" to expand data table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Flake indexing and Slice Height
f_thick = [200,250,300,300,300,     300,300,300,300,    500,   300,300,300,300,300,300,     200]; %Heights of all slices (500 nm is TEM grid)
f_pos = 5000 + cumsum(f_thick) - f_thick/2;
f_pos = f_pos./1e3; %positions of slices as a function of NW length

flake = ["S1F1",string(f_pos(1));
         "S1F2",string(f_pos(2));
         "S1F3",string(f_pos(3));
         "S1F4",string(f_pos(4));
         "S1F5",string(f_pos(5));
         "S2F1",string(f_pos(6));
         "S2F2",string(f_pos(7));
         "S2F3",string(f_pos(8));
         "S2F4",string(f_pos(9));
         "S3F3",string(f_pos(11));
         "S3F4",string(f_pos(12));
         "S3F5",string(f_pos(13));
         "S3F6",string(f_pos(14));
         "S3F7",string(f_pos(15));
         "S3F8",string(f_pos(16));
         "S4F1",string(f_pos(17));
        ];

%% Data Storage (Table Initialisation)
Flake               = cell(16,1); %name of flake
Region_ID           = cell(16,1); %vector for all watershed regions
Height              = cell(16,1);
No_Structures       = cell(16,1);
PL_Spectra          = cell(16,1); %hold PL spectra data points
PL_Fit              = cell(16,1); %hold fitted PL data points
RMS                 = cell(16,1); %RMS of PL fitting
Flagged             = cell(16,1); %whether reion is flagged by unphysical parameters
Region_Size         = cell(16,1); %size of region in pixels
XYZ                 = cell(16,1); %xyz coordinates for most intense pixel per region
Nearest_Region      = cell(16,1); %nearest watershed neighbour
Neighbour_Distance  = cell(16,1); %distance to nearest neighbour
QW_Width            = cell(16,1); %QW Width (interpolated using nextnano results)
Transition_Energy   = cell(16,1);
Wavelength          = cell(16,1);
Fermi_Energy        = cell(16,1);
Urbach_Energy       = cell(16,1);
Sigma               = cell(16,1);
Temperature         = cell(16,1);
Intensity           = cell(16,1); %integrated value of fitted PL spectra curve

QWs_SuperTable      = table(Flake,Region_ID,Height,No_Structures,PL_Spectra,PL_Fit,RMS,Flagged,Region_Size,XYZ,Nearest_Region,...
    Neighbour_Distance,QW_Width,Transition_Energy,Wavelength,Fermi_Energy,Sigma,Temperature,Intensity);

%% Analyse PL Maps
for i = 3 %[11,12,13,14,15,16]
    % 1. Load Data
    disp('1. Loading Data...')

%     load(strcat("..\..\data\PL\",flake(i,1),"_PL_data.mat")); %using locally stored matlab raw data


    settings.wl = h5readatt('..\..\h5 Format\nanoskived_data.h5',strcat('/data/experimental/',flake(i,1),'/pl_map'),'pl_spectrum_x_axis');
    x_range = h5readatt('..\..\h5 Format\nanoskived_data.h5',strcat('/data/experimental/',flake(i,1),'/pl_map'),'map_x_range');
    data = h5read('..\..\h5 Format\nanoskived_data.h5',strcat('/data/experimental/',flake(i,1),'/pl_map'));
    sum_inte = h5read('..\..\h5 Format\nanoskived_data.h5',strcat('/data/experimental/',flake(i,1),'/sum_inte'));

    % 2. Watershed Calculations
    disp('2. Calculating Watershed and Nearest Neighbours...')
    [s,w,regions_found,rgn_size,flatdata,xyz_coords] = PL_watershed(10,1000,settings,data,x_range,sum_inte); %custom script (src)

    % 3. PL Spectra Fitting
    disp('3. Fitting PL Spectra...') 
    regions = 1:regions_found;
    [parameters,PLspectra,PLfit] = PL_LSWfit(settings,w,regions_found,regions,flatdata); %custom script (src)

    % 4. Flagging Unphysical or Weak Data
    disp('4. Flagging Unphysical Parameters')    
    %flag if centre eV not within 750-850 nm
    flag_eV = zeros(regions_found,1);
    flag_eV(regions) = flag_eV(regions) + logical(parameters.Transition_Energy.val(regions) < 1.50);
    flag_eV(regions) = flag_eV(regions) + logical(parameters.Transition_Energy.val(regions) > 1.60);
    flag_eV(regions) = logical(flag_eV);
    eV_PC_reject = 100 * nnz(flag_eV(regions))/length(regions);

    %flag unphysically high temperatures
    flag_temp = zeros(regions_found,1);
    flag_temp(regions) = flag_temp(regions) + logical(parameters.Temperature.val(regions) < 270);
    flag_temp(regions) = flag_temp(regions) + logical(parameters.Temperature.val(regions) > 500);
    flag_temp(regions) = logical(flag_temp);
    temp_PC_reject = 100 * nnz(flag_temp(regions))/length(regions);
    
    %flag if sigma is too large or errorbound < 0 eV
    flag_sigma = zeros(regions_found,1);
    flag_sigma(regions) = flag_sigma(regions) + logical(parameters.Sigma.val(regions) > 0.5);
    flag_sigma(regions) = flag_sigma(regions) + logical(parameters.Sigma.val(regions) < 0);
    flag_sigma(regions) = flag_sigma(regions) + logical(parameters.Sigma.val(regions) - parameters.Sigma.yneg(regions) < 0);
    flag_sigma(regions) = flag_sigma(regions) + logical(parameters.Sigma.val(regions) + parameters.Sigma.ypos(regions) > 0.5);
    flag_sigma(regions) = logical(flag_sigma);
    sigma_PC_reject = 100 * nnz(flag_sigma(regions))/length(regions);

    %flag if QW Width is anonymously large or large errors from LSW fitting
    flag_QW = zeros(regions_found,1);
    flag_QW(regions) = flag_QW(regions) + logical(parameters.QW_Width.val(regions) > 10);
    flag_QW(regions) = flag_QW(regions) + logical(parameters.QW_Width.val(regions) < 3);
    flag_QW(regions) = flag_QW(regions) + logical(parameters.QW_Width.val(regions) - parameters.QW_Width.yneg(regions) < 5);
    flag_QW(regions) = flag_QW(regions) + logical(parameters.QW_Width.val(regions) + parameters.QW_Width.ypos(regions) > 10);
    flag_QW(regions) = logical(flag_QW);
    QW_PC_reject = 100 * nnz(flag_QW(regions))/length(regions);

    %flag anonymously low intensity regions
    flag_int = zeros(regions_found,1);
    flag_int(regions) = flag_int(regions) + isoutlier(parameters.Intensity(regions),'mean');
    flag_int(regions) = logical(flag_int);
    Int_PC_reject = 100 * nnz(flag_int(regions))/length(regions);
    
    %Combine Flags (Pre-NN Calculations)
    o_flag = zeros(regions_found,1); %o(verall) flag filtering regions that aren't evaluated
    o_flag(~regions) = 1;
    p_flag = logical(flag_sigma+flag_eV+flag_temp+flag_QW+flag_int+o_flag);
    
    % 5. Nearest Neighbour
    [n_dist,n_region] = PL_coordinates(regions_found,regions,xyz_coords,p_flag);
    
    %flag Distant Neighbours (to remove dust/scratches etc.)
    flag_NN = zeros(regions_found,1);
    flag_NN(regions) = flag_NN(regions) + logical(n_dist > 5);
    flag_NN = logical(flag_NN);
    NN_PC_reject = 100 * nnz(flag_NN(regions))/length(regions);
    
    %Combine Completed Flagging Stats
    p_flag = logical(flag_sigma+flag_eV+flag_temp+flag_QW+flag_NN+o_flag);
    PCent_Rejected = 100 * nnz(p_flag(regions))/length(regions);
    
    %Summarise Flagging
    Flag_Table = table(p_flag,flag_eV,flag_temp,flag_sigma,flag_NN);
    Flag_Stats = struct('PCent_Rejected',PCent_Rejected,'sigma_Rejected',sigma_PC_reject, ...
        'temp_Rejected',temp_PC_reject,'eV_Rejected',eV_PC_reject,'QW_Rejected',QW_PC_reject, ...
        'Int_Rejected',Int_PC_reject,'NN_Rejected',NN_PC_reject);

    disp(['Regions Flagged for Unphysical TRANSITION ENERGY = ',num2str(Flag_Stats.eV_Rejected),'%'])
    disp(['Regions Flagged for Unphysical SIGMA = ',num2str(Flag_Stats.sigma_Rejected),'%'])
    disp(['Regions Flagged for Unphysical TEMPERATURE = ',num2str(Flag_Stats.temp_Rejected),'%'])
    disp(['Regions Flagged for Unphysical QW WIDTH = ',num2str(Flag_Stats.QW_Rejected),'%'])
    disp(['Regions Flagged for Outlier INTENSITY = ',num2str(Flag_Stats.Int_Rejected),'%'])
    disp(['Regions Flagged for Distant NEAREST NEIGHBOUR = ',num2str(Flag_Stats.NN_Rejected),'%'])
    disp('-----------------------------------------------------------------')
    disp(['Regions Flagged in TOTAL = ',num2str(Flag_Stats.PCent_Rejected),'%'])

    % 6. Organise variables in structures
    disp('6. Sorting Data into Structures...')
    
    rgn_stats = struct('total_regions',regions_found,'rgns_evaluated',regions,'size',rgn_size,'nearest_dist',n_dist,'nearest_region',n_region,'xyz_pos',xyz_coords);
    PLmap = struct('watershed',w,'sum_inte',sum_inte,'x_range',x_range,'y_range',y_range,'region_stats',rgn_stats);

    % 7. Arrange Structure Data in to Tables for each NW regiontemperature.
    disp('7. Arranging into Tables...')
    Region_ID           = (1:regions_found).'; %vector for all watershed regions
    Flagged             = Flag_Table.p_flag; %whether reion is flagged by unphysical parameters
    Region_Size         = rgn_size; %size of region in pixels
    XYZ                 = rgn_stats.xyz_pos; %xyz coordinates for most intense pixel per region
    Nearest_Region      = rgn_stats.nearest_region; %nearest watershed neighbour
    Neighbour_Distance  = rgn_stats.nearest_dist; %distance to nearest neighbour
    QW_Width            = [parameters.QW_Width.val,             parameters.QW_Width.yneg,           parameters.QW_Width.ypos            ];
    Transition_Energy   = [parameters.Transition_Energy.val,    parameters.Transition_Energy.yneg,  parameters.Transition_Energy.ypos   ];
    Wavelength          = [parameters.cow.val,                  parameters.cow.yneg,                parameters.cow.ypos                 ];
    Fermi_Energy        = [parameters.Fermi_Energy.val,         parameters.Fermi_Energy.yneg,       parameters.Fermi_Energy.ypos        ];
    Urbach_Energy       = [parameters.Gamma.val,                parameters.Gamma.yneg,              parameters.Gamma.ypos               ];
    Sigma               = [parameters.Sigma.val,                parameters.Sigma.yneg,              parameters.Sigma.ypos               ];
    Temperature         = [parameters.Temperature.val,          parameters.Temperature.yneg,        parameters.Temperature.ypos         ];
    Intensity           = parameters.Intensity; %integrated value of fitted PL spectra curve
    RMS_Residual        = parameters.rms_r; %residual norm for PL spectra fit
    
    reg_f = Region_ID(~Flagged);

    NWs_Table = table(Region_Size(reg_f),...
        XYZ(reg_f,:),Nearest_Region(reg_f),Neighbour_Distance(reg_f),QW_Width(reg_f,:), ...
        Transition_Energy(reg_f,:),Fermi_Energy(reg_f,:),Urbach_Energy(reg_f,:),Sigma(reg_f,:),Temperature(reg_f,:),...
        Intensity(reg_f),RMS_Residual(reg_f),Wavelength(reg_f,:), ...
        'VariableNames',["Region_Size","XYZ","Nearest_Region","Neighbour_Distance","QW_Width", ...
        "Transition_Energy","Fermi_Energy","Urbach_Energy","Sigma","Temperature","Intensity","RMS_Residual","Wavelength"]);
    PLmap.watershed = uint16(double(ismember(PLmap.watershed,reg_f)).*double(PLmap.watershed));

    % 8. Saving
    disp('8. Saving...')
    save(string(strcat("..\..\results\",flake(i,1),"_PL_Analysed.mat")),'PLmap','NWs_Table','Flag_Table','settings','Flag_Stats')
    disp('9. Saved!')
    
    % 9. Populate Table w/ Data
    QWs_SuperTable.Flake(i)                 = {flake(i,1)};
    QWs_SuperTable.Height(i)                = {str2double(flake(i,2))};
    QWs_SuperTable.No_Structures(i)         = {nnz(~Flagged)};
    QWs_SuperTable.Region_ID(i)             = {Region_ID};
    QWs_SuperTable.PL_Spectra(i)            = {PLspectra};
    QWs_SuperTable.PL_Fit(i)                = {PLfit};
    QWs_SuperTable.RMS(i)                   = {RMS_Residual};       
    QWs_SuperTable.Flagged(i)               = {Flag_Table.p_flag};
    QWs_SuperTable.Region_Size(i)           = {Region_Size};
    QWs_SuperTable.XYZ(i)                   = {XYZ};
    QWs_SuperTable.Nearest_Region(i)        = {Nearest_Region};
    QWs_SuperTable.Neighbour_Distance(i)    = {Neighbour_Distance};
    QWs_SuperTable.QW_Width(i)              = {QW_Width};
    QWs_SuperTable.Transition_Energy(i)     = {Transition_Energy};
    QWs_SuperTable.Wavelength(i)            = {Wavelength};
    QWs_SuperTable.Fermi_Energy(i)          = {Fermi_Energy};
    QWs_SuperTable.Urbach_Energy(i)         = {Urbach_Energy};
    QWs_SuperTable.Sigma(i)                 = {Sigma};
    QWs_SuperTable.Temperature(i)           = {Temperature};
    QWs_SuperTable.Intensity(i)             = {Intensity};
end

%% Save Data in to Expanded Tables
%cell2mat(QWs_SuperTable)
disp('10. Saving All Data...')
QWs_SuperTable.Flake(5) = {string(flake(5,1))}; %Since S1F5 is missing
QWs_SuperTable.Height(5) = {f_pos(5)};
save("..\..\results\QWs_SuperTable_v2.mat","QWs_SuperTable")

disp('11. Saved!')
disp('Completed PL Map Processing.')