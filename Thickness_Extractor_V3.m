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
    
    [obs_results_tot{obs_count}, obs_results_in{obs_count}, obs_results_out{obs_count}, obs_results_cor{obs_count}] = calculate_avg_thickness(obs_folder{obs_count}, spacing, thickness);
    
    % get rid of NaNs
    for i = (1:length(obs_folder{obs_count}))
        obs_results_tot{1,obs_count}(i).avg_thickness_val_right_tot(isnan(obs_results_tot{1,obs_count}(i).avg_thickness_val_right_tot)) = [];
        obs_results_tot{1,obs_count}(i).avg_thickness_val_left_tot(isnan(obs_results_tot{1,obs_count}(i).avg_thickness_val_left_tot)) = [];
        
        obs_results_in{1,obs_count}(i).avg_thickness_val_right_in(isnan(obs_results_in{1,obs_count}(i).avg_thickness_val_right_in)) = [];
        obs_results_in{1,obs_count}(i).avg_thickness_val_left_in(isnan(obs_results_in{1,obs_count}(i).avg_thickness_val_left_in)) = [];
        
        obs_results_out{1,obs_count}(i).avg_thickness_val_right_out(isnan(obs_results_out{1,obs_count}(i).avg_thickness_val_right_out)) = [];
        obs_results_out{1,obs_count}(i).avg_thickness_val_left_out(isnan(obs_results_out{1,obs_count}(i).avg_thickness_val_left_out)) = [];
        
        obs_results_cor{1,obs_count}(i).avg_thickness_val_right_cor(isnan(obs_results_cor{1,obs_count}(i).avg_thickness_val_right_cor)) = [];
        obs_results_cor{1,obs_count}(i).avg_thickness_val_left_cor(isnan(obs_results_cor{1,obs_count}(i).avg_thickness_val_left_cor)) = [];
    end
    
    obs_count = obs_count +1;
    end
end


function [results_tot, results_in, results_out, results_cor] = calculate_avg_thickness(obs_folder, spacing, thickness, obs_count)
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
        seed_px = (left_range + right_range)/2;

        % the segmentation matrix split up into left and right of the seed
        left_of_seed = segmentation(:,1:seed_px);
        % flip the matrix so that the left of the seed can be read left to right
        left_of_seed = flip(left_of_seed,2);
        right_of_seed = segmentation(:,seed_px:end);


        % right of seed
        count2 = 1;
        thick = floor(thickness/2);
        for z = spacing:spacing:size(right_of_seed,2)
            count = 1;
            t_val(count) = right_of_seed(1,z);
            count = count + 1; 
            for y = 1:thick
                t_val(count) = right_of_seed(1,z-y);
                count = count + 1;
                t_val(count) = right_of_seed(1,z+y);
                count = count + 1; 
            end
            t_avg(count2) = mean(t_val);  
            count2 = count2 +1;
        end
        
        
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

            values_left_tot = NaN(1,thickness);
            values_left_in = NaN(1,thickness);
            values_left_out = NaN(1,thickness);
            values_left_cor = NaN(1,thickness);
            
            % for loop to get the values for the thickness sections
            for j = (1:thick_left)
                 values_left_tot(j) = left_of_seed(2,j);   
                 values_left_in(j) = left_of_seed(3,j);
                 values_left_out(j) = left_of_seed(4,j);
                 values_left_cor(j) = left_of_seed(5,j);
            end 
            left_of_seed(:,1:thick_left) = [];
            avg_thickness_val_left_tot(i, count_left) = mean(values_left_tot);
            avg_thickness_val_left_in(i, count_left) = mean(values_left_in);
            avg_thickness_val_left_out(i, count_left) = mean(values_left_out);
            avg_thickness_val_left_cor(i, count_left) = mean(values_left_cor);
            
            count_left = count_left + 1;
        end
        results_tot(i).avg_thickness_val_left_tot = avg_thickness_val_left_tot(i,:);
        results_in(i).avg_thickness_val_left_in = avg_thickness_val_left_in(i,:);
        results_out(i).avg_thickness_val_left_out = avg_thickness_val_left_out(i,:);
        results_cor(i).avg_thickness_val_left_cor = avg_thickness_val_left_cor(i,:);


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

            values_right_tot = NaN(1,thickness);
            values_right_in = NaN(1,thickness);
            values_right_out = NaN(1,thickness);
            values_right_cor = NaN(1,thickness);
            
            % for loop to get the values for the thickness sections
            for j = (1:thick_right)
                 values_right_tot(j) = right_of_seed(2,j);   
                 values_right_in(j) = right_of_seed(3,j);
                 values_right_out(j) = right_of_seed(4,j);
                 values_right_cor(j) = right_of_seed(5,j); 
            end 
            right_of_seed(:,1:thick_right) = [];
            avg_thickness_val_right_tot(i, count_right) = mean(values_right_tot);
            avg_thickness_val_right_in(i, count_right) = mean(values_right_in);
            avg_thickness_val_right_out(i, count_right) = mean(values_right_out);
            avg_thickness_val_right_cor(i, count_right) = mean(values_right_cor);
            count_right = count_right + 1;
        end
        results_tot(i).avg_thickness_val_right_tot = avg_thickness_val_right_tot(i,:);
        results_in(i).avg_thickness_val_right_in = avg_thickness_val_right_in(i,:);
        results_out(i).avg_thickness_val_right_out = avg_thickness_val_right_out(i,:);
        results_cor(i).avg_thickness_val_right_cor = avg_thickness_val_right_cor(i,:);

%         print('Lemme see')

    end
end
















