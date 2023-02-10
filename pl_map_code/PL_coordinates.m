function [n_dist,n_region] = PL_coordinates(regions_found,regions,xyz_coords,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to neighbouring structures for filtering purposes
% (i.e. identifying artefacts that do not correpsond to a heterostructure).
% Calculates co-ordinates of most intense pixel per watershed region
%%% ____________________________ (Inputs) _____________________________ %%%
% regions_found = no. regions found from watershed transform of PL map.
% regions = numbered regions in PL map
% xyz_coords = centre of watershed region
% flag = vector for regions flagged as not containing a heterostructure
%%% ____________________________ (Output) _____________________________ %%%
% n_dist = distance to nearest neighbouring structure
% n_region = identifier (region number) for nearest detected region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_dist = zeros(regions_found,1); %distance to nearest
n_region = zeros(regions_found,1); %nearest region

Region_ID = 1:regions_found;

%rejected regions
rgn_rej = Region_ID(logical(flag));
rgn_rej = rgn_rej(ismember(rgn_rej,regions));


    for i = Region_ID
        %select region coordinates
        selected = [xyz_coords(i,1),xyz_coords(i,2)];
        
        %calculate distance to each region
        reduce = (    xyz_coords(:,1:2) - selected    ).^2;
        dist = sqrt(sum(reduce,2));
        dist(i) = Inf; %make selected region maximum distance (to prevent self-picking)
        dist(rgn_rej) = Inf; %don't include parameter flagged regions
        
        [n_dist(i),n_region(i)] = min(dist(regions));
    end
end