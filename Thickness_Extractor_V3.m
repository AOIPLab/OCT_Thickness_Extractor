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

% select units
list = {'degrees', 'microns', 'pixels'};
[indx_unit,tf1] = listdlg('PromptString','Select Lateral Units', 'SelectionMode', 'single', 'ListString', list);

if tf1 == 1
    if indx_unit ==1
        unit = ' (degrees)';
    elseif indx_unit == 2
        unit = ' (microns)';
    elseif indx_unit == 3
        unit = ' (pixels)';
    end
end

m1 = 'Please enter desired SPACING ';
m2 = 'Please enter the desired WINDOW ';
message1 = strcat(m1, unit);
message2 = strcat(m2, unit);
% specify spacing and thickness
spacing = inputdlg(message1, 'Set Spacing');
spacing = str2double(spacing{1});
thickness = inputdlg(message2, 'Set Sampling Thickness');
thickness = str2double(thickness{1});

% loop to make sure no window overlap
while spacing < (thickness/2)
    warning('Window overlap present with current spacing and window selections');
    contin = questdlg('WARNING: Window overlap with current spacing and window selections Continue with current selection?', ...
      '', ...
      'YES', 'NO', 'NO');
  if strcmpi(contin, 'YES')
    break;
  else
    spacing = inputdlg(message1, 'Set Spacing');
    spacing = str2double(spacing{1});
    thickness = inputdlg(message2, 'Set Sampling Thickness');
    thickness = str2double(thickness{1});
  end
end

% select metrics to have averaged
list = {'Total Retinal Thickness','Inner Retinal Thickness','Outer Retinal Thickness','Corroidal Thickness'};
[indx_metrics,tf2] = listdlg('ListString',list);
    
% Prompt for pop-up asking if you want to add observers
obs_count = 1;
observers = 'YES';
while strcmpi(observers, 'YES')
      observers = questdlg('Add an observer?', ...
      'Add observer', ...
      'YES', 'NO', 'NO');

    if strcmpi(observers, 'NO')
    elseif strcmpi(observers, 'YES')
        % Select observer 2 directory
        obs_fullpath{obs_count} = uigetdir('.','Select directory for observer ');
    
        % Find the files within the chosen directory:
        obs_folder{obs_count} = dir(obs_fullpath{obs_count});
        temp = obs_folder{obs_count};
    
        % Fix the DOS era issue with the dir function (loads in the parent
        % directories '.' and '..')
        obs_folder{obs_count} = temp(~ismember({temp.name}, {'.', '..'}));
        
        % call function
        [obs_results_tot{obs_count}, obs_results_in{obs_count}, obs_results_out{obs_count}, obs_results_cor{obs_count}] = calculate_avg_thickness(obs_folder{obs_count}, spacing, thickness, LUT, indx_unit);
        
        for v = 1:length(obs_folder{1,1})
            % combine left and right values for total
            obs_results_tot{1,obs_count}(v).avg_thickness_total = [obs_results_tot{obs_count}(v).avg_thickness_val_left_tot, obs_results_tot{obs_count}(v).avg_thickness_val_right_tot];
            obs_results_tot{1,obs_count}(v).locations_total = [obs_results_tot{obs_count}(v).locations_left, obs_results_tot{obs_count}(v).locations_right];
            
            % combine left and right values for inner
            obs_results_in{1,obs_count}(v).avg_thickness_total = [obs_results_in{obs_count}(v).avg_thickness_val_left_in, obs_results_in{obs_count}(v).avg_thickness_val_right_in];
            obs_results_in{1,obs_count}(v).locations_total = [obs_results_in{obs_count}(v).locations_left, obs_results_in{obs_count}(v).locations_right];

            % combine left and right values for outer
            obs_results_out{1,obs_count}(v).avg_thickness_total = [obs_results_out{obs_count}(v).avg_thickness_val_left_out, obs_results_out{obs_count}(v).avg_thickness_val_right_out];
            obs_results_out{1,obs_count}(v).locations_total = [obs_results_out{obs_count}(v).locations_left, obs_results_out{obs_count}(v).locations_right];

            % combine left and right values for corroidal
            obs_results_cor{1,obs_count}(v).avg_thickness_total = [obs_results_cor{obs_count}(v).avg_thickness_val_left_cor, obs_results_cor{obs_count}(v).avg_thickness_val_right_cor];
            obs_results_cor{1,obs_count}(v).locations_total = [obs_results_cor{obs_count}(v).locations_left, obs_results_cor{obs_count}(v).locations_right];

        end

        % format output results here
       
        obs_count = obs_count +1;
    end
    
end


function [results_tot, results_in, results_out, results_cor] = calculate_avg_thickness(obs_folder, spacing, thickness, LUT, unit)
    % this for loop goes through each image folder in the grader folder
    % pre-allocate for left
    avg_thickness_val_left = NaN(length(obs_folder), 650);
    % pre-allocate for right
    avg_thickness_val_right = NaN(length(obs_folder), 650);
    
    % initializing results struct
    results_tot = struct([]);
    results_in = struct([]);
    results_out = struct([]);
    results_cor = struct([]);
    
    
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
        
        %add file names to results
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
        
        % find approximate 0 with seed
        % draw seed point
        seed = drawpoint;
        % seed pixel coordinate
        seed_px = round(seed.Position(1));
        
        % find actual 0 from seed reference point
        % look to the left of the seed for the edge
        iteration = 0;
        cor_val1 = segmentation(5,seed_px);
        while cor_val1 == 0
            iteration = iteration+1;
            cor_val1 = segmentation(5,seed_px-iteration);   
        end
        % look to the right of the seed for the edge
        iteration2 = 0;
        cor_val2 = segmentation(5,seed_px);
        while cor_val2 == 0
            iteration2 = iteration2+1;
            cor_val2 = segmentation(5,seed_px+iteration2);
        end
        
        % find the range of the ONH
        left_range = seed_px - iteration+1;
        right_range = seed_px + iteration2-1;
        
        % find actual 0
        % needs to be an integer so round
        seed_px = round((left_range + right_range)/2);

        % the segmentation matrix split up into left and right of the seed
        left_of_seed = segmentation(:,1:seed_px);
        % flip the matrix so that the left of the seed can be read left to right
        left_of_seed = flip(left_of_seed,2);
        right_of_seed = segmentation(:,seed_px:end);
       
        [results_tot, results_in, results_out, results_cor] = calculate(thickness, spacing, left_of_seed, 'left', i, results_tot, results_in, results_out, results_cor, lateral_scale);
        [results_tot, results_in, results_out, results_cor] = calculate(thickness, spacing, right_of_seed, 'right', i, results_tot, results_in, results_out, results_cor, lateral_scale);

    end
end
     

% 
function [results_tot, results_in, results_out, results_cor] = calculate(thickness, spacing, matrix, LorR, i, results_tot, results_in, results_out, results_cor, lateral_scale)

    % convert spacing and window to pixels
    spacing = round(spacing * lateral_scale);
    thickness = round(thickness * lateral_scale);
   
    if mod(thickness,2) % check if odd
        thick_l = floor(thickness/2);
        thick_r = floor(thickness/2);
    else
        thick_l = floor(thickness/2)-1;
        thick_r = floor(thickness/2); % one additional on the right if even
    end

    count2 = 1;
    locations = spacing:spacing:size(matrix,2);
    locations = locations / lateral_scale;
    for z = spacing:spacing:size(matrix,2)
        count = 1;
        values_total(count) = matrix(2,z);
        values_inner(count) = matrix(3,z);
        values_outer(count) = matrix(4,z);
        values_corroidal(count) = matrix(5,z);
        count = count + 1; 
        % go through to get the data on both sides of the sampling
        % point
        for y = 1:thick_l
            values_total(count) = matrix(2,z-y);
            values_inner(count) = matrix(3,z-y);
            values_outer(count) = matrix(4,z-y);
            values_corroidal(count) = matrix(5,z-y);
            count = count + 1;
        end
        for q = 1:thick_r
            % attempt to figure out if not enough data for the whole window
            if (z+q) > size(matrix, 2)
                diff_pos = z+q;
                not_enough = diff_pos - size(matrix, 2);
                incomplete_window = thickness-not_enough;

            else
                values_total(count) = matrix(2,z+q);
                values_inner(count) = matrix(3,z+q);
                values_outer(count) = matrix(4,z+q);
                values_corroidal(count) = matrix(5,z+q);
                count = count + 1; 
                incomplete_window = NaN;
            end
            
        end
        values_total_avg(count2) = mean(values_total);
        values_inner_avg(count2) = mean(values_inner);
        values_outer_avg(count2) = mean(values_outer);
        values_corroidal_avg(count2) = mean(values_corroidal);
        count2 = count2 +1;
    end
    
    % add results to the struct
    if strcmpi(LorR, 'right')
        results_tot(i).avg_thickness_val_right_tot = values_total_avg;
        results_in(i).avg_thickness_val_right_in = values_inner_avg;
        results_out(i).avg_thickness_val_right_out = values_outer_avg;
        results_cor(i).avg_thickness_val_right_cor = values_corroidal_avg;

        results_tot(i).locations_right = locations;
        results_in(i).locations_right = locations;
        results_out(i).locations_right = locations;
        results_cor(i).locations_right = locations;

        results_tot(i).incomplete_window_right = incomplete_window;
        results_in(i).incomplete_window_right = incomplete_window;
        results_out(i).incomplete_window_right = incomplete_window;
        results_cor(i).incomplete_window_right = incomplete_window;

    else
        % flip left results back to read left to right and negate to
        % indicate they are from the left
        results_tot(i).avg_thickness_val_left_tot = flip(values_total_avg,2);
        results_in(i).avg_thickness_val_left_in = flip(values_inner_avg,2);
        results_out(i).avg_thickness_val_left_out = flip(values_outer_avg,2);
        results_cor(i).avg_thickness_val_left_cor = flip(values_corroidal_avg,2);

        results_tot(i).locations_left = -flip(locations);
        results_in(i).locations_left = -flip(locations);
        results_out(i).locations_left = -flip(locations);
        results_cor(i).locations_left = -flip(locations);

        results_tot(i).incomplete_window_left = incomplete_window;
        results_in(i).incomplete_window_left = incomplete_window;
        results_out(i).incomplete_window_left = incomplete_window;
        results_cor(i).incomplete_window_left = incomplete_window;

    end

end
        


















