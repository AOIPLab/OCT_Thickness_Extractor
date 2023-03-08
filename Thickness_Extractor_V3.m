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

%% Loading in LUT file, select lateral units, set spacing and window values, check for overlap, select metricks

% User selects the LUT:
[LUTfile, LUTpath] = uigetfile('.xlsx','Select the LUT');
LUT = readcell(fullfile(LUTpath,LUTfile));

% User selects lateral scale unit
list = {'degrees', 'microns', 'pixels'};
[indx_unit,tf1] = listdlg('PromptString','Select Lateral Units', 'SelectionMode', 'single', 'ListString', list);

% set unit text for spacing and window dailog box
if tf1 == 1
    if indx_unit ==1
        unit = ' (degrees)';
    elseif indx_unit == 2
        unit = ' (microns)';
    elseif indx_unit == 3
        unit = ' (pixels)';
    end
end

% set up messages to be displayed to user to set spacing and window in the desired lateral unit
m1 = 'Please enter desired SPACING ';
m2 = 'Please enter the desired WINDOW ';
message1 = strcat(m1, unit);
message2 = strcat(m2, unit);

% user sets spacing and window
spacing = inputdlg(message1);
spacing = str2double(spacing{1});
window = inputdlg(message2);
window = str2double(window{1});

% spacing must be >= window for no overlap
% if not loop until it is or the user can override and use the overlap
while spacing < (window)
    warning('Window overlap present with current spacing and window selections');
    contin = questdlg('WARNING: Window overlap with current spacing and window selections Continue with current selection?', ...
      '', ...
      'YES', 'NO', 'NO');
  if strcmpi(contin, 'YES')
    break;
  else
    % user sets spacing and window
    spacing = inputdlg(message1);
    spacing = str2double(spacing{1});
    window = inputdlg(message2);
    window = str2double(window{1});
  end
end

% select metrics for output
list = {'Total Retinal Thickness','Inner Retinal Thickness','Outer Retinal Thickness','Corroidal Thickness'};
[indx_metrics,tf2] = listdlg('ListString', list);
    
%% Select observers and call calculation function
% set starting default values
obs_count = 1;
observers = 'YES';

while strcmpi(observers, 'YES')
    % Prompt to ask user to add observers
    observers = questdlg('Add an observer?', ...
    'Add observer', ...
    'YES', 'NO', 'NO');

    if strcmpi(observers, 'NO')
    elseif strcmpi(observers, 'YES')
        % user selects observer directory
        obs_fullpath{obs_count} = uigetdir('.','Select directory for observer ');
    
        % Find the files within the chosen directory:
        obs_folder{obs_count} = dir(obs_fullpath{obs_count});
        temp = obs_folder{obs_count};
    
        % Fix the DOS era issue with the dir function (loads in the parent
        % directories '.' and '..')
        obs_folder{obs_count} = temp(~ismember({temp.name}, {'.', '..'}));
        
        % call 1st function
        [obs_results_tot{obs_count}, obs_results_in{obs_count}, obs_results_out{obs_count}, obs_results_cor{obs_count}] = call_calculate_avg_thickness(obs_folder{obs_count}, spacing, window, LUT, indx_unit);
        
        % data organization
        for v = 1:length(obs_folder{1,1})
            % combine left and right values for total (thickness and locations)
            obs_results_tot{1,obs_count}(v).avg_thickness_total = [obs_results_tot{obs_count}(v).avg_thickness_val_left_tot, obs_results_tot{obs_count}(v).avg_thickness_val_right_tot];
            obs_results_tot{1,obs_count}(v).locations_total = [obs_results_tot{obs_count}(v).locations_left, obs_results_tot{obs_count}(v).locations_right];
            
            % combine left and right values for inner (thickness and locations)
            obs_results_in{1,obs_count}(v).avg_thickness_total = [obs_results_in{obs_count}(v).avg_thickness_val_left_in, obs_results_in{obs_count}(v).avg_thickness_val_right_in];
            obs_results_in{1,obs_count}(v).locations_total = [obs_results_in{obs_count}(v).locations_left, obs_results_in{obs_count}(v).locations_right];

            % combine left and right values for outer (thickness and locations)
            obs_results_out{1,obs_count}(v).avg_thickness_total = [obs_results_out{obs_count}(v).avg_thickness_val_left_out, obs_results_out{obs_count}(v).avg_thickness_val_right_out];
            obs_results_out{1,obs_count}(v).locations_total = [obs_results_out{obs_count}(v).locations_left, obs_results_out{obs_count}(v).locations_right];

            % combine left and right values for corroidal (thickness and locations)
            obs_results_cor{1,obs_count}(v).avg_thickness_total = [obs_results_cor{obs_count}(v).avg_thickness_val_left_cor, obs_results_cor{obs_count}(v).avg_thickness_val_right_cor];
            obs_results_cor{1,obs_count}(v).locations_total = [obs_results_cor{obs_count}(v).locations_left, obs_results_cor{obs_count}(v).locations_right];

        end

        % format output results here
       
        % add to observer count for indexing
        obs_count = obs_count +1;
    end
    
end


%% First funciton - call_calculate_avg_thickness
function [results_tot, results_in, results_out, results_cor] = call_calculate_avg_thickness(obs_folder, spacing, window, LUT, unit)
    
    % pre-allocate for left
    avg_thickness_val_left = NaN(length(obs_folder), 650);
    % pre-allocate for right
    avg_thickness_val_right = NaN(length(obs_folder), 650);
    
    % initializing results structs
    results_tot = struct([]);
    results_in = struct([]);
    results_out = struct([]);
    results_cor = struct([]);
    
    %% get path, load segmentation, find lateral scale, set user seed
    
    % this for loop goes through each image folder in the grader folder
    for i = (1:length(obs_folder))
        
        % get the correct path of where the image and segmentation is
        folder_path = [obs_folder(i).folder '\' obs_folder(i).name];
        image_path = [folder_path '\' obs_folder(i).name '.tif'];
        segmentation_path = [folder_path '\' obs_folder(i).name '_thickness_data.xlsx'];
        
        %load the segmentation
        segmentation = xlsread(segmentation_path, 'thickness_adjusted');  
        
        % open and show the image
        image = imshow(imread(image_path));
        title(obs_folder(i).name, 'Interpreter', 'none');
        
        % add file names to results structs
        results_tot(i).file_name = obs_folder(i).name;
        results_in(i).file_name = obs_folder(i).name;
        results_out(i).file_name = obs_folder(i).name;
        results_cor(i).file_name = obs_folder(i).name;

        % find the index of the lateral scale for the specific scan from
        % the LUT file
        for index = 1:length(LUT)
            trigger = find(strcmp(LUT{index,1},(obs_folder(i).name)));
            if length(trigger) == 1
                break
            end
        end

        % set the lateral scale based on the unit
        if unit == 1 % degrees
            lateral_scale = LUT{index, 2};
        elseif unit == 2 % microns
            lateral_scale = LUT{index, 3};
        elseif unit == 3 % pixels
            lateral_scale = 1;
        end
        
        % find approximate 0 with user selected seed
        % draw seed point
        seed = drawpoint;
        % seed pixel coordinate
        seed_px = round(seed.Position(1));
        
        %%  find actual 0 from seed reference point, split up the segmetation based on left and right of 0, call calculation function
        
        % look to the left of the seed for the edge (end of 0s in corroidal layer)
        iterationl = 0;
        cor_vall = segmentation(5, seed_px);
        while cor_vall == 0
            iterationl = iterationl + 1;
            cor_vall = segmentation(5, seed_px - iterationl);   
        end
        
        % look to the right of the seed for the edge (end of 0s in corroidal layer)
        iterationr = 0;
        cor_valr = segmentation(5, seed_px);
        while cor_valr == 0
            iterationr = iterationr + 1;
            cor_valr = segmentation(5, seed_px+iterationr);
        end
        
        % find the range of the ONH
        left_range = seed_px - iterationl + 1;
        right_range = seed_px + iterationr - 1;
        
        % find actual 0 based on the range - needs to be an integer so round
        seed_px = round((left_range + right_range) / 2);

        % the segmentation matrix split up into left and right of the seed
        left_of_seed = segmentation(:,1:seed_px);
        left_of_seed = flip(left_of_seed, 2); % flip the matrix so that the left of the seed can be read left to right
        right_of_seed = segmentation(:, seed_px:end);
        
        
        % call 2nd function for the left and the right side of 0
        [results_tot, results_in, results_out, results_cor] = calculate_avg_thickness(window, spacing, left_of_seed, 'left', i, results_tot, results_in, results_out, results_cor, lateral_scale);
        [results_tot, results_in, results_out, results_cor] = calculate_avg_thickness(window, spacing, right_of_seed, 'right', i, results_tot, results_in, results_out, results_cor, lateral_scale);

    end
end
     

%% Funciton 2 - Calculate - Does thickness calculations
function [results_tot, results_in, results_out, results_cor] = calculate_avg_thickness(window, spacing, matrix, LorR, i, results_tot, results_in, results_out, results_cor, lateral_scale)

    % convert spacing and window to pixels using the lateral scale - needs to be an integer so round
    spacing = round(spacing * lateral_scale);
    window = round(window * lateral_scale);
   
    if mod(window,2) % check if odd
        window_l = floor(window / 2);
        window_r = floor(window / 2);
    else
        window_l = floor(window / 2) - 1; % subtract 1 so that additional datapoint on the right of the spacing point if even
        window_r = floor(window / 2); 
    end
    
    % get locations stored in an array
    locations = spacing:spacing:size(matrix, 2);
    locations = locations / lateral_scale; % convert locations back to the lateral scale unit to have exact report
    
    spacing_point = 1;
    % for loop to go through all the spacings and calculate thickness for
    % the windows
    for sampling_point = spacing:spacing:size(matrix, 2)
        window_count = 1;
        % get thickness at the sampling point (center of window)
        values_total(window_count) = matrix(2, sampling_point);
        values_inner(window_count) = matrix(3, sampling_point);
        values_outer(window_count) = matrix(4, sampling_point);
        values_corroidal(window_count) = matrix(5, sampling_point);
        window_count = window_count + 1; 
        
        for left_window_p = 1:window_l % go through to get thickness values to the left of the sampling point
            values_total(window_count) = matrix(2, sampling_point - left_window_p);
            values_inner(window_count) = matrix(3, sampling_point - left_window_p);
            values_outer(window_count) = matrix(4, sampling_point - left_window_p);
            values_corroidal(window_count) = matrix(5,sampling_point - left_window_p);
            window_count = window_count + 1;
        end
        
        for right_window_p = 1:window_r % right of sampling point
            % if about to go out of bounds, incomplete window - record the actual size of the window. This would only happen on right side since left matrix gets flipped
            if (sampling_point + right_window_p) > size(matrix, 2)
                pos_diff = sampling_point + right_window_p;
                not_enough = pos_diff - size(matrix, 2);
                incomplete_window = window-not_enough;
                
            else % go through to get thickness values to the right of the sampling point
                values_total(window_count) = matrix(2, sampling_point + right_window_p);
                values_inner(window_count) = matrix(3, sampling_point + right_window_p);
                values_outer(window_count) = matrix(4, sampling_point + right_window_p);
                values_corroidal(window_count) = matrix(5, sampling_point + right_window_p);
                window_count = window_count + 1; 
                incomplete_window = NaN;
            end
            
        end
        % average the values from that sampling window and store
        values_total_avg(spacing_point) = mean(values_total);
        values_inner_avg(spacing_point) = mean(values_inner);
        values_outer_avg(spacing_point) = mean(values_outer);
        values_corroidal_avg(spacing_point) = mean(values_corroidal);
        spacing_point = spacing_point +1;
    end
    
    % add results to the struct
    if strcmpi(LorR, 'right') % if the matrix was to the right of 0
        % store thicknesses
        results_tot(i).avg_thickness_val_right_tot = values_total_avg;
        results_in(i).avg_thickness_val_right_in = values_inner_avg;
        results_out(i).avg_thickness_val_right_out = values_outer_avg;
        results_cor(i).avg_thickness_val_right_cor = values_corroidal_avg;

        % store locations
        results_tot(i).locations_right = locations;
        results_in(i).locations_right = locations;
        results_out(i).locations_right = locations;
        results_cor(i).locations_right = locations;

        % store incomplete window information
        results_tot(i).incomplete_window_right = incomplete_window;
        results_in(i).incomplete_window_right = incomplete_window;
        results_out(i).incomplete_window_right = incomplete_window;
        results_cor(i).incomplete_window_right = incomplete_window;

    else % if the matrix was to the left of 0
        % store thicknesses - flip left results back to read left to right
        results_tot(i).avg_thickness_val_left_tot = flip(values_total_avg,2);
        results_in(i).avg_thickness_val_left_in = flip(values_inner_avg,2);
        results_out(i).avg_thickness_val_left_out = flip(values_outer_avg,2);
        results_cor(i).avg_thickness_val_left_cor = flip(values_corroidal_avg,2);

        % store locations - flip back to read left to right and negate to indicate they are from the left
        results_tot(i).locations_left = -flip(locations);
        results_in(i).locations_left = -flip(locations);
        results_out(i).locations_left = -flip(locations);
        results_cor(i).locations_left = -flip(locations);

        % store incomplete window information
        results_tot(i).incomplete_window_left = incomplete_window;
        results_in(i).incomplete_window_left = incomplete_window;
        results_out(i).incomplete_window_left = incomplete_window;
        results_cor(i).incomplete_window_left = incomplete_window;

    end

end
        
