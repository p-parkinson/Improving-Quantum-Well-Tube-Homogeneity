function [QWs_Energies,matched_wfns] = nextnano_matrix_evaluator(root,QR_name,subno_wfns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate watershed of PL map.
% Calculates co-ordinates of most intense pixel per watershed region
% Following co-ordinates calculation, nearest neighbour distance is
% calculated
%%% ____________________________ (Inputs) _____________________________ %%%
% root = path for simulation output
% QR_name = quantum region name in nextnano quantum{} input file
% subno_wfns = reduce search to this number of wavefunctions
%%% ____________________________ (Output) _____________________________ %%%
% QWs_Energies = 1x3 array of transition energy for the most probable HH-Gamma Transition in QW1,2,3
% matched_wfns = index for correpsonding wavefunctions responsible for transitions in QW1,2,3 in "QWs_Energies"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Import simulation data output and variables
variables_filename = strcat(root,"variables_input.txt"); % variables (e.g. GaAsP Barrier P content)
lastmat_filename = strcat(root,"Structure\last_material_region.origin.dat"); % 2D map of distinctive layered regions
Gamma_filename = strcat(root,"bias_00000\Quantum\probabilities_",QR_name,"_Gamma_.origin.dat"); % Probabilities of Gamma wavefunctions
HH_filename = strcat(root,"bias_00000\Quantum\probabilities_",QR_name,"_HH_.origin.dat"); % Probabilities of HH wavefunctions
Energy_matrix = strcat(root,"bias_00000\Quantum\transition_energies_",QR_name,"_HH_Gamma.origin.dat"); % Energy matrix of HH --> Gamma transitions
Transition_matrix = strcat(root,"bias_00000\Quantum\interband_matrix_elements_",QR_name,"_HH_Gamma.origin.dat"); % HH -> Gamma Transition Probability matrix

%% Last material grid (Mask Creation)
T = readtable(lastmat_filename,'VariableNamingRule','preserve');
minval = min(T.Var2);

grid_length = numel(T.Var2(T.Var2 == minval)); % number of elements in one row/column of square 2D simulation grid
last_mat = reshape(T.Var3,grid_length,grid_length).'; % reshape in to 2D square 

% Mask for QW3 (mask_out), QW2 (mask_mid) and QW1 (mask_inn)
mask_out = last_mat == 5;
mask_mid = last_mat == 7;
mask_inn = last_mat == 9;

%% Load Wavefunctions
variables = readtable(variables_filename,'VariableNamingRule','preserve');
eigens = strcmp(variables.Var1,{'$NumEigenValues'});
index = nonzeros(eigens.*(1:numel(variables.Var1)).');
no_wfns = (variables.Var3(index)); %find no. eigenvalues (wavefunctions)

% Reshape Gamma wavefunctions in to 2D square grid
Gamma = table2array(readtable(Gamma_filename,'VariableNamingRule','preserve'));
rc_elements = Gamma(:,2);
a_length = numel(rc_elements(rc_elements == min(rc_elements)));
Gamma = Gamma(:,3:end);
wfn_G = zeros(a_length,a_length,no_wfns);


HH = table2array(readtable(HH_filename,'VariableNamingRule','preserve'));
HH = HH(:,3:end);
wfn_HH = zeros(a_length,a_length,no_wfns);

% Reshape each Gamma and HH wavefunctions in to 2D square grid (rows,columns)
% in a 3D matrix, where the 3rd dimension corresponds to the index of eigenvalue
for i = 1:no_wfns
    wfn_G(:,:,i) = reshape(Gamma(:,i),a_length,a_length).';
    wfn_HH(:,:,i) = reshape(HH(:,i),a_length,a_length).';
end

%% Bicubic expansion (To match Gamma-HH wavefunctions to QW masks dimensions [i.e. "last_mat"])
adjuster = 0; %floor((numel(mask_out(:,1)) - 2*a_length)/2); %window mask to fit bicubic image
sz = size(last_mat);

%Initialised matrices for expanded Gamma and HH wavefunctions
bc_G = zeros(sz(1),sz(2),no_wfns);
bc_HH = zeros(sz(1),sz(2),no_wfns);

%Expand
for i = 1:no_wfns
    img = imresize(wfn_G(:,:,i),[numel(last_mat(:,1))-adjuster numel(last_mat(1,:))-adjuster],"bicubic"); %weighted avg of 3x3 submatrix
    bc_G(:,:,i) = rescale(img(:,:),0,1);

    img = imresize(wfn_HH(:,:,i),[numel(last_mat(:,1))-adjuster numel(last_mat(1,:))-adjuster],"bicubic"); %weighted avg of 3x3 submatrix
    bc_HH(:,:,i) = rescale(img(:,:),0,1);
end

%% Mask Wavefunctions (Match Gamma-HH modes to QW1,2,3)
% (Uses local function "maxi" which finds the maximum value when overlaying
% 2 arrays (if no overlap, then output will be maxval = 0)

%Arrays for storing in which QW the wavefunctions are present
presence_G = zeros(no_wfns,3); %[in | mid | out]
presence_HH = zeros(no_wfns,3); %[in | mid | out]

for i = 1:no_wfns
    % Gamma
    presence_G(i,1) = maxi(mask_inn,bc_G(:,:,i));
    presence_G(i,2) = maxi(mask_mid,bc_G(:,:,i));
    presence_G(i,3) = maxi(mask_out,bc_G(:,:,i));
    [M,~] = max(presence_G(i,:)); %find which index of QW1-3 is maximum
    presence_G(i,:) = i*double(logical(presence_G(i,:) == M));
    
    % HH
    presence_HH(i,1) = maxi(mask_inn,bc_HH(:,:,i));
    presence_HH(i,2) = maxi(mask_mid,bc_HH(:,:,i));
    presence_HH(i,3) = maxi(mask_out,bc_HH(:,:,i));
    [M,~] = max(presence_HH(i,:)); %find which index of QW1-3 is maximum
    presence_HH(i,:) = i*double(logical(presence_HH(i,:) == M));
end

%% Create masks for the transition matrix for each QW (i.e. exclude wavefunctions > "subno_wfns")
% This also nullifies probabilities for inter-QW transitions (e.g. QW1 -> QW3)
TM_Mask = zeros(100,100,3);
for i = 1:3 %for QW1 to 3
    % Restrict to custom number of wavefunctions (to avoid calculating
    % higher energy states)
    G_modes = presence_G(:,i).*(presence_G(:,i) <= subno_wfns);
    HH_modes = presence_HH(:,i).*(presence_HH(:,i) <= subno_wfns);

    %Create masking matrix
    TM_Mask(:,:,i) = (HH_modes)*(G_modes.') > 0; %[X]columns(HH), [Y]rows(Gamma)
end

%% Append Energies for Maximum Probability Transitions via Overlaying Masked Transitions Matrices to Transition Energy matrix
TX_matrix = readtable(Transition_matrix,'VariableNamingRule','preserve');
TX_matrix = reshape(TX_matrix.Var5,no_wfns,no_wfns);
energies_1D = readtable(Energy_matrix,'VariableNamingRule','preserve');

QWs_Energies = zeros(1,3);
matched_wfns = zeros(3,2);
for i = 1:3 %loop over QW1,2,3
    [Mx,Ix] = max(max(TM_Mask(:,:,i).*TX_matrix,[],2)); %HH
    [My,Iy] = max(max(TM_Mask(:,:,i).*TX_matrix,[],1)); %Gamma
    disp(strcat('Maximum probability: ',num2str(My)))

    %Find Corresponding Energy (Looking at [X,Y] coordinate in Energy Matrix: Var1 = Gamma, Var2 = HH)
    QWs_Energies(i) = abs(nonzeros(energies_1D.Var3.*(energies_1D.Var1 == Ix).*(energies_1D.Var2 == Iy)));

    if Mx == 0 || My == 0
        Ix = 0;
        Iy = 0;
        QWs_Energies(i) = NaN;
        warning(strcat('No Matching Wavefunctions in QW',num2str(i)))
    end
    matched_wfns(i,:) = [Ix,Iy];
end

end % END function


%% Functions
function maxval = maxi(mask,input) %Find maximum value when overlaying 2 arrays (if no overlap, then output will be maxval = 0)
maxval = max(max(mask.*input));
end