function SEM_master = h5_import_sem_master_table(h5_filepath)

flake = h5read(h5_filepath,strcat('/results/SEM_table/','Flake'));
n = height(flake);

SEM_master = table('Size',[n,10],'VariableTypes',...
    ["string","string","double","double","cell", ...
    "double","double","double","double","double"]);
SEM_master.Properties.VariableNames = ["Image","Flake","Height","um_per_px","Centroids",...
    "Diameters","Scale","Nearest Neighbour","Solidity","Eccentricity"];

for vnames = SEM_master.Properties.VariableNames
    var = string(vnames);
%     data = h5read('..\..\h5 Format\nanoskived_results.h5',strcat('/results/PL_table/',var));
%     if ismember(var,["PL_Data","LSW_Fit"])
%         cell_table = cell(n,1);
%         for i = 1:n
%             cell_table(i) = {data(i,:)};
%         end
%         SEM_master.(var) = cell_table;
%     else
        SEM_master.(var) = h5read('..\..\h5 Format\nanoskived_results.h5',strcat('/results/SEM_table/',var));
end

%Convert from string to categorical (h5 doesn't store categorical data)
SEM_master.Flake = categorical(SEM_master.Flake); class(SEM_master.Flake);

end