function master_table = h5_import_pl_master_table(h5_filepath)

flake = h5read(h5_filepath,strcat('/results/PL_table/','Flake'));
n = height(flake);

master_table = table('Size', [n,17],'VariableTypes',...
    ["categorical","cell","cell","double",...
    "table","double","double","double",...
    "logical","double","double","double","double","double" ...
    "double","double","double"]);

master_table.Properties.VariableNames = ["Flake","Height","FlakeID","Size",...
    "XY","PL_Data","LSW_Fit","RMS",...
    "Flag","QW_Width","Transition_Energy","Urbach_Energy","Peak_Emission","Fermi_Energy",...
    "Sigma", "Temperature", "Intensity"];

master_table.PL_Data = cell(n,1);

for vnames = master_table.Properties.VariableNames
    var = string(vnames);
    data = h5read('..\..\h5 Format\nanoskived_results.h5',strcat('/results/PL_table/',var));
    if ismember(var,["PL_Data","LSW_Fit"])
        cell_table = cell(n,1);
        for i = 1:n
            cell_table(i) = {data(i,:)};
        end
        master_table.(var) = cell_table;
    else
        master_table.(var) = h5read('..\..\h5 Format\nanoskived_results.h5',strcat('/results/PL_table/',var));
    end
end

%Convert from string to categorical (h5 doesn't store categorical data)
master_table.Flake = categorical(master_table.Flake); class(master_table.Flake);

end