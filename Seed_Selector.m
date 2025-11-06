% OCT Seed Selector
% Created by: Jenna Grieshop
% Date created: 7/29/2025
%
% Purpose: To allow user to select the seed point for the image so that it
% can be used with the thickness extractor in a separate step
% 
% Requirements: LUT file with the filenames and lateral scale
% information (px/deg, mpp)
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
AxialScale = [];
LateralScale = [];
Seed = [];
LUT_new = table(FileName, AxialScale, LateralScale, Seed);
data = [];

vert_scale_list = {'R2200: 1.6126umpp', 'R2310: 2.4533umpp', 'EnFocus: 5.4076umpp', 'Other (manual entry)'};


% User selects the data directory
data_dir_path = uigetdir('.','Select directory');


% User selects the LUT
try
    [LUTfile, LUTpath] = uigetfile('.csv','Select the Seed LUT, if none - cancel');
    LUT = readtable(fullfile(LUTpath,LUTfile), Delimiter=',');
catch
    new_lut = 1;   
    LUT = table(FileName, AxialScale, LateralScale, Seed);
end

% Get directory of just the .mat files
data_dir_mat = dir(fullfile(data_dir_path, '*.mat'));

% Get list of .mat file names 
for i=1:length(data_dir_mat)
     mat_list{i} = data_dir_mat(i).name;
end
mat_list = string(mat_list);

% if a new LUT is needed, all files are missing, else find what files are
% missing from LUT that exist in the directory
if new_lut == 1 
    missing_from_LUT = mat_list;
else
    lut_list = string(LUT{:,1});
    missing_from_LUT = setdiff(mat_list, lut_list);
end

% Find and add all info to LUT for missing entries
for j=1: length(missing_from_LUT)

    % Get the corresponding image path for the .mat file
    image_path = strrep(missing_from_LUT(j), '.mat', '.tif');

    % open and show the image
    image = imshow(imread(fullfile(data_dir_path, image_path)));
    title(missing_from_LUT(j), 'Interpreter', 'none');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % make the image full screen

    % find approximate 0 with user selected seed
    % draw seed point
    seed = drawpoint;
    % seed pixel coordinate
    seed_px(j) = round(seed.Position(1));
    close

    
    % User to select the vertical scale
    [indx_unit,tf1] = listdlg('PromptString','Select OCT axial(vertical) scale', 'SelectionMode', 'single', 'ListString', vert_scale_list);

    % Set the vertical scale based on user's choice
    if tf1 == 1
        if indx_unit ==1
            vert(j) = 1.6126;
        elseif indx_unit == 2
            vert(j) = 2.4533;
        elseif indx_unit == 3
            vert(j) = 5.4076;
        elseif indx_unit == 4
            m = 'Please enter the axial(vertical) scale in umpp: ';
            temp = inputdlg(m);
            vert(j) = str2double(temp{1});
        end
    end
    
    % Prompt user to enter the lateral scale
    m = 'Please enter the lateral scale in deg/px: ';
    temp2 = inputdlg(m);
    lat(j) = str2double(temp2{1});

  
    data = [data; missing_from_LUT(j), vert(j), lat(j), seed_px(j)];


end

% Add new data to LUT
LUT_new = [LUT; cellstr(data)];
  
% Save LUT
writetable(LUT_new, fullfile(data_dir_path, "Seed_LUT.csv"));  



