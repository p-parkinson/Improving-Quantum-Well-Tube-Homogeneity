%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracts the most probable transition and energy between HH-Gamma
% bandedges at each QW1,2,3 and saves it in a  matrix as a function of QW
% Width (the independent variable of these nextnano simulations). The
% process is carried out for strained systems containing 44%, 47% and 50% P
% content within the GaAsP barriers and an unstrained system with P = 47%
% for comparison.
% This is based on the script 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

loops = 3:0.5:12; %QW Widths simulated
P_comps = [0.44,0.47,0.5]; %P contents for strained simulations

QWs_Energy_loop = zeros(numel(loops),3,numel(P_comps)+1); % +1 to append unstrained values
QR_name = "quantum_region";

%% Start Script
u_root = strcat("D:\OUTPUT\(Unstrained) 2D GaAsP_GaAs Radial QW Heterostructures_QWThickness_");
for k = 1:numel(P_comps)+1
    
    for j = 1:numel(loops)
        
        if k < numel(P_comps)+1 %strained
            s_root = strcat("D:\OUTPUT\(Aruni) 2D GaAsP_GaAs Radial QW Heterostructures_GaAsPAlloy_Barrier_",num2str(P_comps(k)),"_QWThickness_");
            if ismember(loops(j),3:12)
                f_root = strcat(s_root,num2str(loops(j)),".0\");
            else
                f_root = strcat(s_root,num2str(loops(j)),"\");
            end
            QWs_Energy_loop(j,:,k) = nextnano_matrix_evaluator(f_root,QR_name,0,8);
        else %unstrained
            if ismember(loops(j),3:12)
                f_root = strcat(u_root,num2str(loops(j)),".0\");
            else
                f_root = strcat(u_root,num2str(loops(j)),"\");
            end
            QWs_Energy_loop(j,:,k+1) = nextnano_matrix_evaluator(f_root,QR_name,0,50);
        end
    end

end

Strained_Energies = QWs_Energy_loop(:,:,1:3);
Unstrained_Energies = QWs_Energy_loop(:,:,end);

disp('DONE')