% OCT Thickness Extractor
% Created by: Jenna Grieshop
% Date created: 3/9/2023
%
% Purpose: Interactive script to calculate average OCT thicknesses from one or more 
% graders/observers at specific sampling points (spacing) with a specific
% bin (window) sizes from the center of the optic nerve head
% 
% Requirements: LUT file with the file/folder names and lateral scale
% information (px/deg, mpp), grader/observer folder(s)
%
% Inputs: prompt for LUT file, prompt for unit selection, prompt for
% spacing and window size, prompt to select observer/grader folder, prompt
% to select "seed" on image of approximate 0 (center of the ONH)
%
% Outputs: An excel file for each grader will be generated in the location
% of the LUT file.
%
% Note: If the window size will overlap with the selections made
% it will give a warning and another chance to set the sampling and window
% size.


clear all
close all
clc

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

% select metrics for output - commented out for now and hard coded to just
% use total and choroidal - if want to use this in the future would need to
% finish the code by adding in which ones to use in the output data below
% list = {'Total Retinal Thickness','Inner Retinal Thickness','Outer Retinal Thickness','Choroidal Thickness'};
% [indx_metrics,tf2] = listdlg('ListString', list);

% header for output file
header = {'File Name', 'Location', 'Total Retinal Thickness', 'Choroidal Thickness', 'Inner Thickness', 'Outer Thickness' 'Bin Size'};
    
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
        data_compiled_1 = [];
        % user selects observer directory
        obs_fullpath{obs_count} = uigetdir('.','Select directory for observer ');
    
        % Find the files within the chosen directory:
        obs_folder{obs_count} = dir(obs_fullpath{obs_count});
        temp = obs_folder{obs_count};
    
        % Fix the DOS era issue with the dir function (loads in the parent
        % directories '.' and '..')
        obs_folder{obs_count} = temp(~ismember({temp.name}, {'.', '..'}));

        obs_name = split(obs_fullpath, '\');
        obs_name = obs_name(end);
        output_fname = strcat(obs_name{1,1}, '_OCT_Thickness_Results_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')));        
        
        % call 1st function
        [obs_results_tot{obs_count}, obs_results_in{obs_count}, obs_results_out{obs_count}, obs_results_chor{obs_count}] = call_calculate_avg_thickness(obs_folder{obs_count}, spacing, window, LUT, indx_unit);
        
        % data organization
        for v = 1:length(obs_folder{1,obs_count})
            
            % combine 0s for total thickness
            obs_results_tot{1,obs_count}(v).zero_thickness_total = [obs_results_tot{obs_count}(v).zero_vals_left_tot, obs_results_tot{obs_count}(v).zero_vals_right_tot]; % combine the left and right 0 data
            l_r_0_mean = mean(obs_results_tot{1,obs_count}(v).zero_thickness_total); % get the mean of the left and right 0s thicknesses
            obs_results_tot{obs_count}(v).avg_thickness_val_left_tot(end) = l_r_0_mean; % add mean thickness value to the end of left
            obs_results_tot{obs_count}(v).avg_thickness_val_right_tot(1) = []; % remove first thickness element from right
            obs_results_tot{obs_count}(v).size_of_bin_left(end) = obs_results_tot{obs_count}(v).size_of_bin_left(end) + obs_results_tot{obs_count}(v).size_of_bin_right(1); % add the 0 bin sizes
            obs_results_tot{obs_count}(v).size_of_bin_right(1) = [];% remove first size of bin element from right
            obs_results_tot{obs_count}(v).locations_right(1) = []; % remove first location element from right
            
            % combine 0s for choroidal thickness
            obs_results_chor{1,obs_count}(v).zero_thickness_total = [obs_results_chor{obs_count}(v).zero_vals_left_chor, obs_results_chor{obs_count}(v).zero_vals_right_chor]; % combine the left and right 0 data
            l_r_0_mean_chor = mean(obs_results_chor{1,obs_count}(v).zero_thickness_total); % get the mean of the left and right 0s thicknesses
            obs_results_chor{obs_count}(v).avg_thickness_val_left_chor(end) = l_r_0_mean_chor; % add mean thickness value to the end of left
            obs_results_chor{obs_count}(v).avg_thickness_val_right_chor(1) = []; % remove first thickness element from right
            obs_results_chor{obs_count}(v).size_of_bin_left(end) = obs_results_chor{obs_count}(v).size_of_bin_left(end) + obs_results_chor{obs_count}(v).size_of_bin_right(1); % add the 0 bin sizes
            obs_results_chor{obs_count}(v).size_of_bin_right(1) = [];% remove first size of bin element from right
            obs_results_chor{obs_count}(v).locations_right(1) = []; % remove first location element from right

            % combine 0s for inner thickness
            obs_results_in{1,obs_count}(v).zero_thickness_total = [obs_results_in{obs_count}(v).zero_vals_left_in, obs_results_in{obs_count}(v).zero_vals_right_in]; % combine the left and right 0 data
            l_r_0_mean_in = mean(obs_results_in{1,obs_count}(v).zero_thickness_total); % get the mean of the left and right 0s thicknesses
            obs_results_in{obs_count}(v).avg_thickness_val_left_in(end) = l_r_0_mean_in; % add mean thickness value to the end of left
            obs_results_in{obs_count}(v).avg_thickness_val_right_in(1) = []; % remove first thickness element from right
            obs_results_in{obs_count}(v).size_of_bin_left(end) = obs_results_in{obs_count}(v).size_of_bin_left(end) + obs_results_in{obs_count}(v).size_of_bin_right(1); % add the 0 bin sizes
            obs_results_in{obs_count}(v).size_of_bin_right(1) = [];% remove first size of bin element from right
            obs_results_in{obs_count}(v).locations_right(1) = []; % remove first location element from right

            % combine 0s for outer thickness
            obs_results_out{1,obs_count}(v).zero_thickness_total = [obs_results_out{obs_count}(v).zero_vals_left_out, obs_results_out{obs_count}(v).zero_vals_right_out]; % combine the left and right 0 data
            l_r_0_mean_out = mean(obs_results_out{1,obs_count}(v).zero_thickness_total); % get the mean of the left and right 0s thicknesses
            obs_results_out{obs_count}(v).avg_thickness_val_left_out(end) = l_r_0_mean_out; % add mean thickness value to the end of left
            obs_results_out{obs_count}(v).avg_thickness_val_right_out(1) = []; % remove first thickness element from right
            obs_results_out{obs_count}(v).size_of_bin_left(end) = obs_results_out{obs_count}(v).size_of_bin_left(end) + obs_results_out{obs_count}(v).size_of_bin_right(1); % add the 0 bin sizes
            obs_results_out{obs_count}(v).size_of_bin_right(1) = [];% remove first size of bin element from right
            obs_results_out{obs_count}(v).locations_right(1) = []; % remove first location element from right

            % combine left and right values for total (thickness, locations, bin sizes)
            obs_results_tot{1,obs_count}(v).avg_thickness_total = [obs_results_tot{obs_count}(v).avg_thickness_val_left_tot, obs_results_tot{obs_count}(v).avg_thickness_val_right_tot];
            obs_results_tot{1,obs_count}(v).locations_total = [obs_results_tot{obs_count}(v).locations_left, obs_results_tot{obs_count}(v).locations_right];
            obs_results_tot{1,obs_count}(v).size_of_bin_total = [obs_results_tot{obs_count}(v).size_of_bin_left, obs_results_tot{obs_count}(v).size_of_bin_right];

            % combine left and right values for inner (thickness, locations, bin sizes)
            obs_results_in{1,obs_count}(v).avg_thickness_total = [obs_results_in{obs_count}(v).avg_thickness_val_left_in, obs_results_in{obs_count}(v).avg_thickness_val_right_in];
            obs_results_in{1,obs_count}(v).locations_total = [obs_results_in{obs_count}(v).locations_left, obs_results_in{obs_count}(v).locations_right];
            obs_results_in{1,obs_count}(v).size_of_bin_total = [obs_results_in{obs_count}(v).size_of_bin_left, obs_results_in{obs_count}(v).size_of_bin_right];

            % combine left and right values for outer (thickness, locations, bin sizes)
            obs_results_out{1,obs_count}(v).avg_thickness_total = [obs_results_out{obs_count}(v).avg_thickness_val_left_out, obs_results_out{obs_count}(v).avg_thickness_val_right_out];
            obs_results_out{1,obs_count}(v).locations_total = [obs_results_out{obs_count}(v).locations_left, obs_results_out{obs_count}(v).locations_right];
            obs_results_out{1,obs_count}(v).size_of_bin_total = [obs_results_out{obs_count}(v).size_of_bin_left, obs_results_out{obs_count}(v).size_of_bin_right];

            % combine left and right values for choroidal (thickness, locations, bin sizes)
            obs_results_chor{1,obs_count}(v).avg_thickness_total = [obs_results_chor{obs_count}(v).avg_thickness_val_left_chor, obs_results_chor{obs_count}(v).avg_thickness_val_right_chor];
            obs_results_chor{1,obs_count}(v).locations_total = [obs_results_chor{obs_count}(v).locations_left, obs_results_chor{obs_count}(v).locations_right];
            obs_results_chor{1,obs_count}(v).size_of_bin_total = [obs_results_chor{obs_count}(v).size_of_bin_left, obs_results_chor{obs_count}(v).size_of_bin_right];

            % formatting output for each grader

            % get the name of the scan in a cell array with the necessary
            % ammount
            file_name = cell(length(obs_results_tot{1,obs_count}(v).locations_total),1);
            for p = 1:length(file_name)
                file_name{p,1} = obs_results_tot{obs_count}(v).file_name;
            end
            
            % Combine the data with the file names and previous data from the same observer
            if length(data_compiled_1) == 0
                data_compiled_1 = [file_name num2cell(obs_results_tot{1,obs_count}(v).locations_total') num2cell(obs_results_tot{1,obs_count}(v).avg_thickness_total') num2cell(obs_results_chor{1,obs_count}(v).avg_thickness_total') num2cell(obs_results_in{1,obs_count}(v).avg_thickness_total') num2cell(obs_results_out{1,obs_count}(v).avg_thickness_total') num2cell(obs_results_tot{1,obs_count}(v).size_of_bin_total')]; 
            else
                data_compiled_1 = [data_compiled_1; file_name num2cell(obs_results_tot{1,obs_count}(v).locations_total') num2cell(obs_results_tot{1,obs_count}(v).avg_thickness_total') num2cell(obs_results_chor{1,obs_count}(v).avg_thickness_total') num2cell(obs_results_in{1,obs_count}(v).avg_thickness_total') num2cell(obs_results_out{1,obs_count}(v).avg_thickness_total') num2cell(obs_results_tot{1,obs_count}(v).size_of_bin_total')];
            end

        end
        % combine headers with the data
        data_compiled_2 = [header; data_compiled_1];
        % write output file for the observer
        xlswrite(fullfile(LUTpath,output_fname), data_compiled_2);

        % add to observer count for indexing
        obs_count = obs_count +1;
    end
    
end
close all

%% First funciton - call_calculate_avg_thickness
function [results_tot, results_in, results_out, results_chor] = call_calculate_avg_thickness(obs_folder, spacing, window, LUT, unit)
    
    % pre-allocate for left
    avg_thickness_val_left = NaN(length(obs_folder), 650);
    % pre-allocate for right
    avg_thickness_val_right = NaN(length(obs_folder), 650);
    
    % initializing results structs
    results_tot = struct([]);
    results_in = struct([]);
    results_out = struct([]);
    results_chor = struct([]);
    
    %% get path, load segmentation, find lateral scale, set user seed
    
    % this for loop goes through each image folder in the grader folder
    for i = (1:length(obs_folder))
        
        % get the correct path of where the image and segmentation is
        folder_path = [obs_folder(i).folder '\' obs_folder(i).name];
        image_path = [folder_path '\' obs_folder(i).name '.tif'];
        segmentation_path = [folder_path '\' obs_folder(i).name '_thickness_data.xlsx'];
        
        %load the segmentation
        segmentation = xlsread(segmentation_path, 'thickness_raw');  
        
        % open and show the image
        image = imshow(imread(image_path));
        title(obs_folder(i).name, 'Interpreter', 'none');
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % make the image full screen
        
        % add file names to results structs
        results_tot(i).file_name = obs_folder(i).name;
        results_in(i).file_name = obs_folder(i).name;
        results_out(i).file_name = obs_folder(i).name;
        results_chor(i).file_name = obs_folder(i).name;

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
        chor_vall = segmentation(5, seed_px);
        while chor_vall == 0
            iterationl = iterationl + 1;
            chor_vall = segmentation(5, seed_px - iterationl);   
        end
        
        % look to the right of the seed for the edge (end of 0s in corroidal layer)
        iterationr = 0;
        chor_valr = segmentation(5, seed_px);
        while chor_valr == 0
            iterationr = iterationr + 1;
            chor_valr = segmentation(5, seed_px+iterationr);
        end
        
        % find the range of the ONH
        left_range = seed_px - iterationl + 1;
        right_range = seed_px + iterationr - 1;
        
        % find actual 0 based on the range - needs to be an integer so round
        turd = (left_range + right_range) / 2;
        seed_px = round((left_range + right_range) / 2);

        % the segmentation matrix split up into left and right of the seed
        left_of_seed = segmentation(:,1:seed_px);
        left_of_seed = flip(left_of_seed, 2); % flip the matrix so that the left of the seed can be read left to right
        right_of_seed = segmentation(:, seed_px+1:end);
        
        
        % call 2nd function for the left and the right side of 0
        [results_tot, results_in, results_out, results_chor] = calculate_avg_thickness(window, spacing, left_of_seed, 'left', i, results_tot, results_in, results_out, results_chor, lateral_scale);
        [results_tot, results_in, results_out, results_chor] = calculate_avg_thickness(window, spacing, right_of_seed, 'right', i, results_tot, results_in, results_out, results_chor, lateral_scale);

    end
end
     
%% Funciton 2 - Calculate - Does thickness calculations
function [results_tot, results_in, results_out, results_cor] = calculate_avg_thickness(window, spacing, matrix, LorR, i, results_tot, results_in, results_out, results_cor, lateral_scale)
 
    % initialize the bin_size value so that if the bin size remains 0 the code won't try to add non-existant data to the struct
    bin_size = 0;

    size_of_matrix_in_unit = size(matrix, 2) * (1/lateral_scale); % get the size of the matrix in desired unit
    spacing_starting_point = 0; % have the starting point at 0
    sampling_points = spacing_starting_point:spacing:size_of_matrix_in_unit; % get locations of all the sampling points in desired unit
    if strcmpi(LorR, 'left')
        px_to_unit = 0:(1/lateral_scale):50; % pixel location converted to unit value to max degrees possible
    else
        px_to_unit = (1/lateral_scale):(1/lateral_scale):50; % pixel location converted to unit value to max degrees possible
    end
    px_to_unit = px_to_unit(1:size(matrix,2)); % pixel location converted to unit value at each index corresponding to the data matrix

    window_count = 1;
    first = 1;
    first_store = 1;
    for point = sampling_points

        if first
            % get the min and max for the bin range
            bin_range_min = point;
            bin_range_max = point + (window/2);
            first = 0;
        else
            % get the min and max for the bin range
            bin_range_min = point - (window/2);
            bin_range_max = point + (window/2);
        end
        
        px_index = 1;
        count = 1;
        indices = [];
        % go through each unit converted pixel location to determine which indeces in the matrix are in the window
        for location = px_to_unit
            if location >= bin_range_min % location must be greater than or equal to the bin range minimum
                if location < bin_range_max % location must be less than the bin range maximum
                    indices(count) = px_index; % store the index of the location that fits in the bin range
                    count = count +1;
                end
            end
            px_index = px_index +1;
        end
        
        % find the index range for values that fit in the window
        index_max = max(indices);
        index_min = min(indices);
        % get the ammount of data points that fit in the window
        bin_size(window_count) = length(indices);
        % average the thickness of the matrix indices that fit in the window and store the values for each window
        values_total(window_count) = mean(matrix(2, (index_min:index_max)));
        values_inner(window_count) = mean(matrix(3, (index_min:index_max)));
        values_outer(window_count) = mean(matrix(4, (index_min:index_max)));
        values_choroidal(window_count) = mean(matrix(5, (index_min:index_max)));

        % store the raw values from the 0 bin to average later
        if first_store
            values_total_zero = matrix(2, (index_min:index_max));
            values_inner_zero = matrix(3, (index_min:index_max));
            values_outer_zero = matrix(4, (index_min:index_max));
            values_choroidal_zero = matrix(5, (index_min:index_max));
            first_store = 0;
        end

        window_count = window_count +1;

    end
    
    % add results to the struct
    if strcmpi(LorR, 'right') % if the matrix was to the right of 0
        if bin_size == 0 % if bin size is zero there is no data to store
        else
            % store thicknesses
            results_tot(i).zero_vals_right_tot = values_total_zero;
            results_tot(i).avg_thickness_val_right_tot = values_total;
            results_in(i).zero_vals_right_in = values_inner_zero;
            results_in(i).avg_thickness_val_right_in = values_inner;
            results_out(i).zero_vals_right_out = values_outer_zero;
            results_out(i).avg_thickness_val_right_out = values_outer;
            results_cor(i).zero_vals_right_chor = values_choroidal_zero;
            results_cor(i).avg_thickness_val_right_chor = values_choroidal;
    
            % store locations
            results_tot(i).locations_right = sampling_points;
            results_in(i).locations_right = sampling_points;
            results_out(i).locations_right = sampling_points;
            results_cor(i).locations_right = sampling_points;
    
            % store incomplete window information
            results_tot(i).size_of_bin_right = bin_size;
            results_in(i).size_of_bin_right = bin_size;
            results_out(i).size_of_bin_right = bin_size;
            results_cor(i).size_of_bin_right = bin_size;
        end

    else % if the matrix was to the left of 0
        if bin_size == 0 % if bin size is zero there is no data to store
        else
        % store thicknesses - flip left results back to read left to right
        results_tot(i).zero_vals_left_tot = flip(values_total_zero);
        results_tot(i).avg_thickness_val_left_tot = flip(values_total,2);
        results_in(i).zero_vals_left_in = flip(values_inner_zero);
        results_in(i).avg_thickness_val_left_in = flip(values_inner,2);
        results_out(i).zero_vals_left_out = flip(values_outer_zero);
        results_out(i).avg_thickness_val_left_out = flip(values_outer,2);
        results_cor(i).zero_vals_left_chor = flip(values_choroidal_zero);
        results_cor(i).avg_thickness_val_left_chor = flip(values_choroidal,2);

        % store locations - flip back to read left to right and negate to indicate they are from the left
        results_tot(i).locations_left = -flip(sampling_points);
        results_in(i).locations_left = -flip(sampling_points);
        results_out(i).locations_left = -flip(sampling_points);
        results_cor(i).locations_left = -flip(sampling_points);

        % store incomplete window information
        results_tot(i).size_of_bin_left = flip(bin_size);
        results_in(i).size_of_bin_left = flip(bin_size);
        results_out(i).size_of_bin_left = flip(bin_size);
        results_cor(i).size_of_bin_left = flip(bin_size);
        end

    end

end
        
