% OCT Seed Selector
% Created by: Jenna Grieshop
% Date created: 7/29/2025
%
% Purpose: To allow user to select the seed point for the image so that it
% can be used with the thickness extractor in a separate step
% 
% Requirements: LUT file with the filenames and lateral scale
% information (px/deg, mpp). If you do not have this, select cancel when
% prompted and a new LUT will be generated from this script.
%
% Inputs: prompt for LUT file, prompt for unit selection, prompt
% to select "seed" on image of approximate 0 (center of the ONH)
%
% Outputs: An excel file with the seed point will be generated or if one
% exists in the folder the information will be added to the LUT file
%


clear all
close all
clc
new_lut = 0;
FileName = strings(0);
AxialScaleUMPP = [];
LateralScalePPD = [];
LateralScaleUMPP = [];
Seed = [];
PixelSizeX = [];
PixelSizeY = [];
NominalSizeX = [];
NominalSizeY = [];
AL = [];
LUT_new = table(FileName, AxialScaleUMPP, LateralScalePPD, LateralScaleUMPP, Seed, PixelSizeX, PixelSizeY, NominalSizeY);
data = [];
dtype = 1;

vert_scale_list = {'R2200: 1.6126umpp', 'R2310: 2.4533umpp', 'EnFocus: 5.4076umpp', 'Other (manual entry)'};


% User selects the data directory
data_dir_path = uigetdir('.','Select directory');

% Get directory of just the .mat files
data_dir = dir(fullfile(data_dir_path, '*.mat'));
if isempty(data_dir)
    data_dir = dir(fullfile(data_dir_path, '*.csv'));
    containsOutput = contains({data_dir.name}, 'LUT');
    data_dir = data_dir(~containsOutput);
    dtype = 2;
    data = table(FileName, AxialScaleUMPP, LateralScalePPD, LateralScaleUMPP, Seed, AL, PixelSizeX, PixelSizeY, NominalSizeX, NominalSizeY);
else
    data = table(FileName, AxialScaleUMPP, LateralScalePPD, LateralScaleUMPP, Seed);
end


% User selects the LUT
try
    [LUTfile, LUTpath] = uigetfile('.csv','Select the Seed LUT, if none - cancel');
    LUT = readtable(fullfile(LUTpath,LUTfile), Delimiter=',');
catch
    new_lut = 1;  
    if dtype == 1
        LUT = table(FileName, AxialScaleUMPP, LateralScalePPD, LateralScaleUMPP, Seed);
    else
        LUT = table(FileName, AxialScaleUMPP, LateralScalePPD,  LateralScaleUMPP, Seed, AL, PixelSizeX, PixelSizeY, NominalSizeX, NominalSizeY);
    end
end

% Get list of data file names 
for i=1:length(data_dir)
     data_list{i} = data_dir(i).name;
end
data_list = string(data_list);

% if a new LUT is needed, all files are missing, else find what files are
% missing from LUT that exist in the directory
if new_lut == 1 
    missing_from_LUT = data_list;
    missing_seed_name = [];
else
    lut_list = string(LUT{:,1});
    missing_from_LUT = setdiff(data_list, lut_list);
    nan_indices = isnan(LUT.Seed);
    missing_seed = LUT(nan_indices, :);
    missing_seed_name = missing_seed.FileName;
end

% Find and add all info to LUT for missing entries
for j=1: length(missing_from_LUT)

    if dtype == 1
        % Get the corresponding image path for the .mat file
        image_path = strrep(missing_from_LUT(j), '.mat', '.tif');
    else
        % Get the corresponding image path for the .csv file
        image_path = strrep(missing_from_LUT(j), '.csv', '.tif');
    end

    % open and show the image
    image = imshow(imread(fullfile(data_dir_path, image_path)));
    [px_y(j), px_x(j), chn] = size(image.CData);
    title(missing_from_LUT(j), 'Interpreter', 'none');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % make the image full screen

    % find approximate 0 with user selected seed
    % draw seed point
    seed = drawpoint;
    % seed pixel coordinate
    seed_px(j) = round(seed.Position(1));
    close

    
    % User to select the vertical scale
    [indx_unit,tf1] = listdlg('PromptString','Select OCT axial (vertical) scale', 'SelectionMode', 'single', 'ListString', vert_scale_list);

    % Set the vertical scale based on user's choice
    if tf1 == 1
        if indx_unit ==1
            vert(j) = 1.6126;
        elseif indx_unit == 2
            vert(j) = 2.4533;
        elseif indx_unit == 3
            vert(j) = 5.4076;
        elseif indx_unit == 4
            m = 'Please enter the axial (vertical) scale in microns per pixel for ' + missing_from_LUT(j) + ': ';
            temp = inputdlg(m);
            vert(j) = str2double(temp{1});
        end
    end
    
    
    % Prompt user to enter the lateral scale ppd
    m = 'Please enter the lateral scale in pixel per degree (ppd) for ' + missing_from_LUT(j) + ': ';
    temp2 = inputdlg(m);
    lat_ppd(j) = str2double(temp2{1});

    % Prompt user to enter the lateral scale umpp
    m = 'Please enter the lateral scale in microns per pixel (umpp) for ' + missing_from_LUT(j) + ': ';
    temp3 = inputdlg(m);
    lat_umpp(j) = str2double(temp3{1});
    nominal_y(j) = vert(j) * px_y(j);
    nominal_x(j) = lat_umpp(j) * px_x(j);

    if dtype == 2
        
        %Prompt user to enter the Axial Length
        m = 'Please enter the axial length for '+ missing_from_LUT(j) + ': ';
        temp4 = inputdlg(m);
        AL(j) = str2double(temp4{1});

        data = [data; {missing_from_LUT(j), vert(j), lat_ppd(j), lat_umpp(j) seed_px(j), AL(j), px_x(j), px_y(j), nominal_x(j), nominal_y(j)}];
    else
        data = [data; {missing_from_LUT(j), vert(j), lat_ppd(j), lat_umpp(j), seed_px(j)}];


    end

end

%% fill in rows with missing seed points only
for k = 1:length(missing_seed_name)
     if dtype == 1
        % Get the corresponding image path for the .mat file
        image_path = strrep(missing_from_LUT(k), '.mat', '.tif');
    else
        % Get the corresponding image path for the .csv file
        image_path = strrep(missing_from_LUT(k), '.csv', '.tif');
    end

    % open and show the image
    image = imshow(imread(fullfile(data_dir_path, image_path)));
    [px_y(k), px_x(k), chn] = size(image.CData);
    title(missing_from_LUT(k), 'Interpreter', 'none');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % make the image full screen

    % find approximate 0 with user selected seed
    % draw seed point
    seed = drawpoint;
    % seed pixel coordinate
    seed_px(k) = round(seed.Position(1));
    close

    index = strcmp(LUT.FileName, missing_seed.FileName);
    LUT{index, 'Seed'} = seed_px(k);

end

if ~isempty(missing_from_LUT)
    % Add new data to LUT
    LUT_new = [LUT; data];
else
    LUT_new = LUT;
end
  
% Save LUT
writetable(LUT_new, fullfile(data_dir_path, strcat("Seed_LUT_", string(datetime('now','TimeZone','local','Format','yyyyMMdd')), ".csv")));  



