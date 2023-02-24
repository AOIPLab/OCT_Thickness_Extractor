% Thickness Extractor
% Jenna's Take

clear all
close all
clc

%% Load in the lib folder 
mfilefullpath = matlab.desktop.editor.getActiveFilename; 
mfileshortpath = erase(mfilefullpath, 'Thickness_Extractor_V1.m'); 
lib_path = append(mfileshortpath, 'lib');
addpath(genpath(lib_path));

%% Loading in LUT file, OCT image folder, Observer segmentation spreasheets 

% User selects the LUT:
[LUTfile, LUTpath] = uigetfile('.xlsx','Select the LUT');
LUT = readcell(fullfile(LUTpath,LUTfile));

% User selects the desired directory containing all images:
img_fullpath = uigetdir('.','Select directory containing images');

% Find the files within the chosen directory:
img_folder = dir(img_fullpath);

% Fix the DOS era issue with the dir function (loads in the parent
% directories '.' and '..')
img_folder = img_folder(~ismember({img_folder.name}, {'.', '..'}));

% Pre- allocate cell array for all images
OCT_images = cell(length(img_folder),1);

% Select observer 1 directory
obs1_fullpath = uigetdir('.','Select directory for observer 1 ');

% Find the files within the chosen directory:
obs1_folder = dir(obs1_fullpath);

% Fix the DOS era issue with the dir function (loads in the parent
% directories '.' and '..')
obs1_folder = obs1_folder(~ismember({obs1_folder.name}, {'.', '..'}));

% initializing results struct
obs1_results = struct([]);


% select units
list = {'degrees', 'um', 'pixels'};
[indx_unit,tf1] = listdlg('ListString',list);

% specify spacing and thickness
spacing = inputdlg('Please enter the desired SPACING', 'Set Spacing');
spacing = round(str2double(spacing{1}));
thickness = inputdlg('Please enter the desired SAMPLING THICKNESS', 'Set Thickness');
thickness = round(str2double(thickness{1}));

% select metrics to have averaged
list = {'Total Retinal Thickness','Inner Retinal Thickness','Outer Retinal Thickness','Corroidal Thickness'};
[indx_metrics,tf2] = listdlg('ListString',list);

obs1_results = calculate_avg_thickness(obs1_folder, obs1_results, spacing, thickness);
    

% probs put this in a function as well and just have a global counter going
% on to keep track of what number grader we are on
% Prompt for pop-up asking if you have more than one observer
obs_count = 1;
while observers ~= 'NO'
    observers = questdlg('Add an observer?', ...
      'Add observer', ...
      'YES', 'NO', 'NO');

    if strcmpi(observers, 'NO')
    elseif strcmpi(observers, 'YES')
    % Select observer 2 directory
    obs_fullpath(obs_count) = uigetdir('.','Select directory for observer 2 ');

    % Find the files within the chosen directory:
    obs_folder(obs_count) = dir(obs_fullpath(obs_count));
    temp = obs_folder(obs_count);

    % Fix the DOS era issue with the dir function (loads in the parent
    % directories '.' and '..')
    obs_folder(obs_count) = temp(~ismember({obs2_folder.name}, {'.', '..'}));

    % initializing results struct
    obs_results(obs_count) = struct([]);

    obs_count = obs_count +1;
    end
end


function obs_results = calculate_avg_thickness(obs_folder, obs_results, spacing, thickness)
    % possibly put this in a function so that when each new grader is added
    % this gets called again
    % this for loop goes through each image folder in the grader folder
    % pre-allocate for left
    avg_thickness_val_left = NaN(length(obs_folder), 650);
    % pre-allocate for right
    avg_thickness_val_right = NaN(length(obs_folder), 650);
    for i = (1:length(obs_folder))
        % get the correct path of where the image and segmentation is
        folder_path = [obs_folder(i).folder '\' obs_folder(i).name];
        image_path = [folder_path '\' obs_folder(i).name '.tif'];
        segmentation_path = [folder_path '\' obs_folder(i).name '_thickness_data.xlsx'];
        %load the segmentation
        segmentation = xlsread(segmentation_path, 'thickness_adjusted');    
        % open and show the image
        image = imshow(imread(image_path));
        % draw seed point
        seed = drawpoint;
        % seed pixel coordinate
        seed_px = round(seed.Position(1));

        % the segmentation matrix split up into left and right of the seed
        left_of_seed = segmentation(:,1:seed_px-1);
        % flip the matrix so that the left of the seed can be read left to
        % right
        left_of_seed = flip(left_of_seed,2);
        right_of_seed = segmentation(:,seed_px+1:end);


        % should put this as a function and have left/right_of_seed as an
        % input
        % left of the seed
        count_left = 1;
        while size(left_of_seed,2) > 0

            if spacing > 0
                if spacing > size(left_of_seed,2)
                    left_of_seed(:,1:size(left_of_seed,2)) = [];
                else
                    left_of_seed(:,1:spacing) = [];
                end
            end


            % if less than the thickness is left, just use what is left
            if thickness < size(left_of_seed,2)
                thick_left = thickness;
            else
                thick_left = size(left_of_seed,2);
            end

            values_left = NaN(1,thickness);
            % for loop to get the values for the thickness sections
            for j = (1:thick_left)
                 values_left(j) = left_of_seed(2,j); % would need to store all values and then choose later which to be reported   
            end 
            left_of_seed(:,1:thick_left) = [];
            avg_thickness_val_left(i, count_left) = mean(values_left);
            count_left = count_left + 1;
        end
        obs_results(i).avg_thickness_val_left = avg_thickness_val_left(i,:);


        % right of the seed
        count_right = 1;
        while size(right_of_seed,2) > 0

            if spacing > 0
                if spacing > size(right_of_seed,2)
                    right_of_seed(:,1:size(right_of_seed,2)) = [];
                else
                    right_of_seed(:,1:spacing) = [];
                end
            end


            % if less than the thickness is left, just use what is left
            if thickness < size(right_of_seed,2)
                thick_right = thickness;
            else
                thick_right = size(right_of_seed,2);
            end

            values_right = NaN(1,thickness);
            % for loop to get the values for the thickness sections
            for j = (1:thick_right)
                 values_right(j) = right_of_seed(2,j); % would need to store all values and then choose later which to be reported   
            end 
            right_of_seed(:,1:thick_right) = [];
            avg_thickness_val_right(i, count_right) = mean(values_right);
            count_right = count_right + 1;
        end
        obs_results(i).avg_thickness_val_right = avg_thickness_val_right(i,:);

        print('Lemme see')

    end
end
















