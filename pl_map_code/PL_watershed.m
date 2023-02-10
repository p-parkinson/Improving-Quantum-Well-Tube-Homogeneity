 function [s,w,regions_found,rgn_size,flatdata,xyz_coords] = PL_watershed(px_low,px_hi,settings,data,x_range,sum_inte)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate watershed of PL map.
% Calculates co-ordinates of most intense pixel per watershed region
% Following co-ordinates calculation, nearest neighbour distance is
% calculated
%%% ____________________________ (Inputs) _____________________________ %%%
% px_low = lower bound for no. pixels per watershed region
% px_hi = upper bound for no. pixels per watershed region
% settings = import settings from measurement
% data = PL spectra
% coe = peak emission of each pixel
% x_range = co-ordinates of square PL map (-50um to 50 um, 334 steps)
%%% ____________________________ (Output) _____________________________ %%%
% s = size of PL map
% w = watershed transform
% regions_found = number of regions found from watershed
% rgn_size = size of each watershed region in no. pixels
% flatdata = PL Spectrum data (reshaped to correct dimension)
% xyz_coords = centre of watershed region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts=optimset('Display','off'); %suppress fitting messages
%% Initialisation
x = settings.wl; %setup x-axis (wavelength, nm)
xr = 1e3.*x_range;
n = @(x)((x-min(x(:)))/(max(x(:))-min(x(:))));
g = @(x,p)(p(1).*exp(-(x-p(2)).^2/(2*p(3).^2)));

% Reshape spectra
flatdata = reshape(data,[],1024);
s = size(flatdata);
    %% Prepare Fitting for Watershed (Based on PL Spectra)
    fitp = zeros(s(1),4);

    %% Fit Spectra in Preparation for Watershed Transform
    %  (to find peak intensity and remove low intensity spectra)
    for i = 1:s(1)
        % Get data
        y = double(flatdata(i,:));
        % Sort from big to small
        [s,o] = sort(y,'descend');
        % Peak as median of 10 largest values
        mh = median(s(1:10));
        % Position as median of position of 10 largest values
        mp = median(x(o(1:10)));
        % Exclude data with low intensity peak
        if mh < 20
            fitp(i,4) = 0;
            continue;
        end
        % Function to minimize
        tomin = @(p)(sum((y-mean(y(1:200))-g(x,p)).^2));
        % Do fit
        [p,r] = fminsearch(tomin,[mh, mp, 600],opts);
        % Save best fit values
        fitp(i,1:3) = p;
        % Save residual
        fitp(i,4) = r;
    end
    % Cut out bad fits
    d_mask = fitp(:,4)>1e7;
    fitp(d_mask,1:3) = 0;

    %% Watershed 
    % Calculate watershed (to separate wires based on intensity)
    s = size(sum_inte);
    w = watershed(-reshape(fitp(:,1),s(1),s(2)),8);
    regions_found = max(w(:));
    disp_r = regions_found;
    disp(['Detected ',num2str(disp_r),' Regions.'])
    
    %%% Remove regions that are small, update to chronological numbering.
    j = 1; %counter
    for i = 1:regions_found
        %Examine isolated regions
        single_nw = logical(w==i);
        sz = nnz(single_nw);
        if sz > px_low && sz < px_hi
            %Update watershed
            w(w==i) = j;
            %update region number
            j = j + 1;
        else
            w(w==i) = 0;
        end
    end
    
    regions_found = max(w(:)); %update watershed regions

%% Find Co-ordinates of Each Region
xyz_coords = zeros(regions_found,3);
rgn_size = zeros(regions_found,1); %size of region in no. pixels
    
    for q = 1:regions_found
        rgn_size(q) = nnz(logical(w==q));
        single_nw = logical(w==q).*sum_inte';
        
        %Find coordinates at peak intensity
        [diam_max,diam_indx] = max(max(single_nw));
        [~,diam_indy] = max(max(single_nw.'));
        xyz_coords(q,:) = [xr(diam_indx) xr(diam_indy) diam_max];
    end
    
end