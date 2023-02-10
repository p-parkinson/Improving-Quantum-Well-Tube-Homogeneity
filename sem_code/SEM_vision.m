%% Analyse SEM data to obtain diameters etc.

% Flake metadata
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

%% Settings and Folders
S1_tifs = dir("../../data/SEM/S1/*.tif");
S2_tifs = dir("../../data/SEM/S2/*.tif");
S3_tifs = dir("../../data/SEM/S3/*.tif");
S4_tifs = dir("../../data/SEM/S4/*.tif");

sub = ["S1","S2","S3","S4"];
%% Start Data Parsing
for j = 1:4 %loop for 1-4 substrates
    clear data
    if j == 1
        tifs = S1_tifs;
    elseif j == 2
        tifs = S2_tifs;
    elseif j == 3
        tifs = S3_tifs;
    elseif j == 4
        tifs = S4_tifs;
    else
        error('ERROR at assigning tifs')
    end
    f_path = {tifs.folder}.'; %List all paths, length = number of images per substrate
    i_name = {tifs.name}.';
    
    %%% DATA CELL FOR STORING IMAGING STATS ETC.
    data = cell(numel(f_path),15);

    %Assignment of ambiguous image names to SXFY format and assign height
    if j == 1
        data(1:5,2) = {"S1F1"};
        data([6:12,42,43],2) = {"S1F2"};
        data(7:26,2) = {"S1F3"};
        data([27:41,44:49],2) = {"S1F4"};
        data(50,2) = {"sample"};
        
        data(1:5,15) = {num2str(f_pos(1))};
        data([6:12,42,43],15) = {num2str(f_pos(2))};
        data(7:26,15) = {num2str(f_pos(3))};
        data([27:41,44:49],15) = {num2str(f_pos(4))};
        data(50,2) = {"NA"};
    elseif j == 2
        data(1:8,2) = {"S2F1"};
        data(9:17,2) = {"S2F2"};
        data(18:31,2) = {"S2F3"};
        data(32:45,2) = {"S2F4"};
        data(46,2) = {"sample"};
        
        data(1:8,15) = {num2str(f_pos(6))};
        data(9:17,15) = {num2str(f_pos(7))};
        data(18:31,15) = {num2str(f_pos(8))};
        data(32:45,15) = {num2str(f_pos(9))};
        data(46,2) = {"NA"};
    elseif j == 3
        data(2:8,2) = {"S3F3"};
        data(9:15,2) = {"S3F4"};
        data(16:19,2) = {"S3F5"};
        data(20:24,2) = {"S3F6"};
        data(25:28,2) = {"S3F7"};
        data(29:34,2) = {"S3F8"};
        data([1,35,36],2) = {"sample"};
        
        data(2:8,15) = {num2str(f_pos(11))};
        data(9:15,15) = {num2str(f_pos(12))};
        data(16:19,15) = {num2str(f_pos(13))};
        data(20:24,15) = {num2str(f_pos(14))};
        data(25:28,15) = {num2str(f_pos(15))};
        data(29:34,15) = {num2str(f_pos(16))};
        data([1,35,36],15) = {"sample"};
    elseif j == 4
        data(:,2) = {"S4F1"};
        
        data(:,15) = {num2str(f_pos(17))};
    else
        disp('ERROR at assigning flakes')
    end

    %% Text Recognition for Metadata
    bar_pxlen = zeros(numel(f_path),1); %pixel length of scale bar
    sb_length = zeros(numel(f_path),1); %scale bar length in microns
    um_per_px = zeros(numel(f_path),1); %resolution
    for i =1:numel(f_path)
        i_path = strcat(f_path(i),"\",i_name(i));
        im_SEM = imread(i_path,"tif");
        im_meta = imbinarize(imcrop(im_SEM,[1,888,1022,55]));
    
            % Find px length of scale bar
            imbar = imcrop(im_meta,[493 1 531 27]);
            stats = regionprops(imbar,'Area','BoundingBox','Centroid','MajorAxisLength');
            MAL = [stats.MajorAxisLength]; %filter by major axis length (MAL)
            BB = reshape([stats.BoundingBox],[4,numel(stats)]);
            MAL = (MAL<400 & MAL>30 & BB(2,:)==9.5);
            stats = stats(MAL==1);
            b_box = [stats.BoundingBox];
            if size(b_box,2) ~= 8
                error('Misdetection of scale bar regions')
            end
            bar_pxlen(i) = abs((b_box(5)+b_box(7))-b_box(1)); %find range of bounding box to get bar length
    
            % Scale bar numeric length extract
            ocr = ocr(im_meta,[714 1 88 26]);
            d = strjoin(regexp(string(ocr.Words(1)),'([0-9]|[1-9]0|[1-9]00)','match')); %extract number
            sb_length(i) = str2double(regexprep(d, '\s+', '')); %remove deadspace
            clear ocr
            
            data(i,3) = {sb_length(i)};
            %Save SEM image without metadata
            SEM_Image = imcrop(im_SEM,[0,0,1024,883]);
            %imwrite(SEM_Image,strcat("SEM_Cropped/",string(data(i,2)),"/",string(data(i,2)),"_",string(data(i,3)),"um_(",num2str(i),").png"))
    end
    um_per_px(:) = sb_length(:)./bar_pxlen(:);
    sb_table = table(sb_length,um_per_px);
    
    %% Regionprops (Load 2 images and find XY coordinates)
    min_r = 150; %min radius [nm] of object
    max_r = 2000; %max radius [nm] of object
    for i = 1:numel(f_path)
        i_path = strcat(f_path(i),"\",i_name(i));
        immat = imread(i_path); %read image in to matlab
        imgray = imcrop(immat,[0 0 1024 885]); %crop out metadata
        immat = imbinarize(imgray); %binarize image for regionprops
    
        %Store data in table
        data(i,1) = {strcat(f_path(i),"\",i_name(i))}; %file path
        data(i,3) = {sb_length(i)}; % scale bar length (e.g. 10 um)
        data(i,4) = {um_per_px(i)}; % resolution
        
            % Analyse images with MATLAB's regionprops function
            stats = regionprops('table',immat,'Centroid',...
                'MajorAxisLength','MinorAxisLength','Eccentricity','Solidity','Orientation');
            data(i,5) = {stats.Centroid}; % xy coordinate of object
            data(i,6) = {mean([stats.MajorAxisLength stats.MinorAxisLength],2)*um_per_px(i)}; % diameters [um]
            data(i,7) = {mean([stats.MajorAxisLength stats.MinorAxisLength],2)*1e3*um_per_px(i)}; % diameters [nm]
            data(i,8) = {stats.MajorAxisLength*1e3*um_per_px(i) < max_r & stats.MinorAxisLength*1e3*um_per_px(i) > min_r}; % flag accepted values (1) if within 150-2000nm diameter
            %Diameter analysis
            diams = cell2mat(data(i,7));
            flags = logical(cell2mat(data(i,8)));
            data(i,9) = {mean(diams(flags))}; %mean diameter
            data(i,10) = {std(diams(flags))}; %std diameter
            data(i,13) = {stats.Eccentricity}; %eccentricity
            data(i,14) = {stats.Solidity}; %solidity
            
            regions_found = numel(data{i,5}(:,1));
            regions = 1:regions_found;
            regions_found = regions(data{i,8});
    end
    
    %% Nearest Neighbours (NNs)
    for i = 1:numel(f_path)
        regions_found = numel(data{i,5}(:,1));
        n_rgns = regions_found;
        regions = 1:regions_found;
        regions_found = regions(data{i,8});
        
        n_dist = zeros(n_rgns,1); %distance to nearest
        n_region = zeros(n_rgns,1); %nearest region
        dist = zeros(n_rgns,1);
        xyz_coords = data{i,5};
    
        if nnz(regions) < numel(n_rgns) %if flagged data is present
            xyz_coords(~regions,:) = [Inf,inf]; %make flagged data Inf (x,y) so it isn't picked up as a NN
        end
        
        for k = regions_found
            %select region coordinates
            selected = [xyz_coords(k,1),xyz_coords(k,2)];
    
            %calculate distance to each region
            reduce = (    xyz_coords(:,1:2) - selected    ).^2;
            dist = sqrt(sum(reduce,2));
            dist(k) = Inf; %make selected region maximum distance (to prevent self-picking)
    
            [n_dist(k),n_region(k)] = min(dist(:));
        end
        data(i,11) = {n_dist(:).*data{i,4}};
        data(i,12) = {n_region(:)};
    end
    
    %% Saving
    disp('Saving...')

    %shorten file path (so universal for other systems)
    %THIS IS UNIQUE TO FULL PATH!!split(path,"\");
    for i = 1:height(data(:,1))
        [path,name,ext] = fileparts(string(data(i,1)));
        path = split(path,"\");
        data(i,1) = {strcat("..\data\SEM\",path(end),"\",name,ext)};
    end
    if j == 1
        S1_SEM = struct('filename',data(:,1),'flake',data(:,2),'Scale',data(:,3),'um_per_px',data(:,4),...
        'centroids',data(:,5),'diameters_px',data(:,6),'diameters_nm',data(:,7),'flag',data(:,8),...
        'Mean_Diameter',data(:,9),'STD_Diameter',data(:,10),'NN_Dist',data(:,11),'NN_RegionID',data(:,12),...
        'Eccentricity',data(:,13),'Solidity',data(:,14),'FlakeHeight',data(:,15));
    elseif j == 2
        S2_SEM = struct('filename',data(:,1),'flake',data(:,2),'Scale',data(:,3),'um_per_px',data(:,4),...
        'centroids',data(:,5),'diameters_px',data(:,6),'diameters_nm',data(:,7),'flag',data(:,8),...
        'Mean_Diameter',data(:,9),'STD_Diameter',data(:,10),'NN_Dist',data(:,11),'NN_RegionID',data(:,12),...
        'Eccentricity',data(:,13),'Solidity',data(:,14),'FlakeHeight',data(:,15));
    elseif j == 3
        S3_SEM = struct('filename',data(:,1),'flake',data(:,2),'Scale',data(:,3),'um_per_px',data(:,4),...
        'centroids',data(:,5),'diameters_px',data(:,6),'diameters_nm',data(:,7),'flag',data(:,8),...
        'Mean_Diameter',data(:,9),'STD_Diameter',data(:,10),'NN_Dist',data(:,11),'NN_RegionID',data(:,12),...
        'Eccentricity',data(:,13),'Solidity',data(:,14),'FlakeHeight',data(:,15));
    elseif j == 4
        S4_SEM = struct('filename',data(:,1),'flake',data(:,2),'Scale',data(:,3),'um_per_px',data(:,4),...
        'centroids',data(:,5),'diameters_px',data(:,6),'diameters_nm',data(:,7),'flag',data(:,8),...
        'Mean_Diameter',data(:,9),'STD_Diameter',data(:,10),'NN_Dist',data(:,11),'NN_RegionID',data(:,12),...
        'Eccentricity',data(:,13),'Solidity',data(:,14),'FlakeHeight',data(:,15));
    else
        error('ERROR at saving')
    end
end
SEM_Data = [S1_SEM;S2_SEM;S3_SEM;S4_SEM];
SEM_Data = struct2table(SEM_Data);
save(string(strcat('..\..\results\SEM_Data.mat')),'SEM_Data')

%% Create SEM Master Table
no_structures = [];
for r = 1:height(SEM_Data)
    no_structures = [no_structures; nnz(cell2mat(SEM_Data{r,"flag"}))];
end

n=sum(no_structures);
SEM_master = table('Size',[n,10],'VariableTypes',...
    ["string","string","double","double","cell", ...
    "double","double","double","double","double"]);
SEM_master.Properties.VariableNames = ["Image","Flake","Height","um_per_px","Centroids",...
    "Diameters","Scale","Nearest Neighbour","Solidity","Eccentricity"];

centroids_table = table('Size',[n,2],'VariableTypes',["double","double"]);
centroids_table.Properties.VariableNames = ["X","Y"];

L = 1; %counter
for r = 1:height(SEM_Data)
    mask = cell2mat(SEM_Data{r,"flag"});

    dx = nnz(mask);
    
    SEM_master{L:L+dx-1,"Image"} = SEM_Data{r,"filename"};
    SEM_master{L:L+dx-1,"Flake"} = SEM_Data{r,"flake"};
    SEM_master{L:L+dx-1,"Scale"} = SEM_Data{r,"Scale"};
    SEM_master{L:L+dx-1,"um_per_px"} = SEM_Data{r,"um_per_px"};
    centroids_table{L:L+dx-1,:} = SEM_Data{r,'centroids'}{1}(mask,:);
    SEM_master{L:L+dx-1,"Diameters"} = SEM_Data{r,"diameters_nm"}{1}(mask);
    SEM_master{L:L+dx-1,"Nearest Neighbour"} = SEM_Data{r,"NN_Dist"}{1}(mask);
    SEM_master{L:L+dx-1,"Height"} = str2double(SEM_Data{r,"FlakeHeight"});
    SEM_master{L:L+dx-1,"Solidity"} = SEM_Data{r,"Solidity"}{1}(mask);
    SEM_master{L:L+dx-1,"Eccentricity"} = SEM_Data{r,"Eccentricity"}{1}(mask);

    L = L+dx;
end

SEM_master.Flake = categorical(SEM_master.Flake);
SEM_master.Centroids = table2array(centroids_table);

%Filter out unwanted data
SEM_master = SEM_master((SEM_master.Flake ~= "sample"),:);
SEM_master = SEM_master((SEM_master.Flake ~= "NA"),:);
SEM_master = SEM_master((SEM_master.Scale < 40),:); %don't analyse images with low resolution
SEM_master = SEM_master((SEM_master.Solidity > 0.9),:);

save(string(strcat('..\..\results\SEM_master.mat')),'SEM_master')