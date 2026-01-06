% OCT Thickness Extractor
% Created by: Jenna Grieshop
% Date created: 3/9/2023
%
% Purpose: Interactive script to calculate average OCT thicknesses from
% .mat files generated in doctrap. This uses specific sampling points (spacing) 
% with specific bin (window) sizes from the center of the optic nerve head
% (ONH location determined from Seed_Selector.m script)
% 
% Requirements: LUT file generated from the Seed_Selector.m script
%
% Inputs: prompt for LUT file, prompt for unit selection, prompt for
% spacing and window size, prompt for data folder
%
% Outputs: An excel file with results will be generated in the location
% of the LUT file. Rows are as follows: eccentricity locations, number of
% datapoints contributing to the average, average thickness between all
% data at that eccentricity, average thickness for each file.
%
% Note: If the window size will overlap with the selections made
% it will give a warning and another chance to set the sampling and window
% size.


clear all
close all
clc

basePath = which('Thickness_Extractor_new.m');

[basePath ] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); % Add our support library to the path.

%% Loading in LUT file, select lateral units, set spacing and window values, check for overlap, select metricks

% User selects the LUT:
[LUTfile, LUTpath] = uigetfile('.csv','Select the LUT');
LUT = readcell(fullfile(LUTpath,LUTfile));

% User selects lateral scale unit
list = {'degrees', 'microns', 'pixels'};
[indx_unit,tf1] = listdlg('PromptString','Select Lateral Units', 'SelectionMode', 'single', 'ListString', list);

% set unit text for spacing and window dailog box
if tf1 == 1
    if indx_unit ==1
        unit = 'degrees';
    elseif indx_unit == 2
        unit = 'microns';
    elseif indx_unit == 3
        unit = 'pixels';
    end
end

% set up messages to be displayed to user to set spacing and window in the desired lateral unit
m1 = 'Please enter desired SPACING in ';
m2 = 'Please enter the desired WINDOW in ';
message1 = strcat(m1, unit);
message2 = strcat(m2, unit);

% user sets spacing and window
spacing_str = inputdlg(message1);
spacing = str2double(spacing_str{1});
window_str = inputdlg(message2);
window = str2double(window_str{1});

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
    spacing_str = inputdlg(message1);
    spacing = str2double(spacing_str{1});
    window_str = inputdlg(message2);
    window = str2double(window_str{1});
  end
end



%% User selects data directory

% user selects data directory
data_dir_path = uigetdir('.','Select directory of data ');



% Get directory of just the .mat files (used for DOCTRAP version)
data_dir = dir(fullfile(data_dir_path, '*.mat'));
if isempty(data_dir)

    % Get directory of just the .csv files (used for generic version)
    data_dir = dir(fullfile(data_dir_path, '*.csv'));
    dtype = 2; % for csv general data
else
    dtype = 1; % for .mat doctrap data
end

% here is where we need to check the file type in the directory and then
% set a flag for either using doctrap or not - JG generalize DONE 1/5/26

% Get list of data file names 
for i=1:length(data_dir)
     data_list{i} = data_dir(i).name;
end
data_list = string(data_list);

%% Call functions to calculate thickness and organize

% call function to calculate avg thicknesses
[results, results_tot, results_in, results_out, results_chor] = call_calculate_avg_thickness(dtype, data_dir, spacing, window, LUT, indx_unit);

% data organization - call to combine_results function
results_tot = combine_results(results_tot, data_dir);
results_in = combine_results(results_in, data_dir);
results_out = combine_results(results_out, data_dir);
results_chor = combine_results(results_chor, data_dir);

%% Get max and min locations between the datasets

% initialize min and max location
min_loc = 0;
max_loc = -10;

% find overall min and max locations for the datasets
for v = 1:length(data_dir)
    if min_loc > min(results_tot(v).locations_left)
                min_loc = min(results_tot(v).locations_left);
    end
    if max_loc < max(results_tot(v).locations_right)
        max_loc = max(results_tot(v).locations_right);
    end
end

%% Sort data for output

% Get the step size betweeen the bin centers
all_bin_locs = min_loc:spacing:max_loc;

% Initialize results tables
results_values_tot = zeros(length(all_bin_locs),length(data_dir));
results_sum_tot = zeros(length(all_bin_locs),1);
results_contributed_tot = zeros(length(all_bin_locs),1);

results_values_in = zeros(length(all_bin_locs),length(data_dir));
results_sum_in = zeros(length(all_bin_locs),1);
results_contributed_in = zeros(length(all_bin_locs),1);

results_values_out = zeros(length(all_bin_locs),length(data_dir));
results_sum_out = zeros(length(all_bin_locs),1);
results_contributed_out = zeros(length(all_bin_locs),1);

results_values_chor = zeros(length(all_bin_locs),length(data_dir));
results_sum_chor = zeros(length(all_bin_locs),1);
results_contributed_chor = zeros(length(all_bin_locs),1);

% Find which subjects have data at each possible location
for z = 1:length(all_bin_locs')
    clear results_list;
    for j = 1:length(data_dir)
        max_size_of_bin = max(results_tot(j).size_of_bin_total); 
        name_list{j} =  results_tot(j).file_name;
        for k = 1:length(results_tot(j).locations_total)
            if results_tot(j).locations_total(k) == all_bin_locs(z)
                if results_tot(j).size_of_bin_total(k) >= (max_size_of_bin - 1) % make sure to only use full bins
 
                    % total thickness
                    results_values_tot(z,j) = results_tot(j).avg_thickness_total(k);
                    results_sum_tot(z) = results_sum_tot(z) + results_tot(j).avg_thickness_total(k);
                    results_contributed_tot(z) = results_contributed_tot(z) + 1;
                    
                    % inner thickness
                    results_values_in(z,j) = results_in(j).avg_thickness_total(k);
                    results_sum_in(z) = results_sum_in(z) + results_in(j).avg_thickness_total(k);
                    results_contributed_in(z) = results_contributed_in(z) + 1;
                    
                    % outer thickness
                    results_values_out(z,j) = results_out(j).avg_thickness_total(k);
                    results_sum_out(z) = results_sum_out(z) + results_out(j).avg_thickness_total(k);
                    results_contributed_out(z) = results_contributed_out(z) + 1;
                    
                    % choroidal thickness
                    results_values_chor(z,j) = results_chor(j).avg_thickness_total(k);
                    results_sum_chor(z) = results_sum_chor(z) + results_chor(j).avg_thickness_total(k);
                    results_contributed_chor(z) = results_contributed_chor(z) + 1;
    
                    break
                end
            end
        end    
    end
end

%% format data for output and save

header = ['Eccentricity', 'Datapoints Contributed', 'Average', name_list]';
data_compiled_tot = [num2cell(all_bin_locs); num2cell(results_contributed_tot)'; num2cell((results_sum_tot./results_contributed_tot))';  num2cell(results_values_tot)'];
data_compiled_total = [header, data_compiled_tot];

data_compiled_in = [num2cell(all_bin_locs); num2cell(results_contributed_in)'; num2cell((results_sum_in./results_contributed_in))';  num2cell(results_values_in)'];
data_compiled_inner = [header, data_compiled_in];

data_compiled_out = [num2cell(all_bin_locs); num2cell(results_contributed_out)'; num2cell((results_sum_out./results_contributed_out))';  num2cell(results_values_out)'];
data_compiled_outer = [header, data_compiled_out];

data_compiled_chor = [num2cell(all_bin_locs); num2cell(results_contributed_chor)'; num2cell((results_sum_chor./results_contributed_chor))';  num2cell(results_values_chor)'];
data_compiled_choroidal = [header, data_compiled_chor];


% replace . with p for better file naming convention
spacing_str = strrep(spacing_str, '.', 'p');
window_str = strrep(window_str, '.', 'p');

% write output file for the observer
output_fname = strcat('OCT_Thickness_Results_spacing_', spacing_str{1}, '_window_', window_str{1}, '_', unit, '_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.xlsx');        
writecell(data_compiled_total, fullfile(LUTpath,output_fname), 'Sheet', 'Total Thickness');
writecell(data_compiled_inner, fullfile(LUTpath,output_fname), 'Sheet', 'Inner Thickness');
writecell(data_compiled_outer, fullfile(LUTpath,output_fname), 'Sheet', 'Outer Thickness');
writecell(data_compiled_choroidal, fullfile(LUTpath,output_fname), 'Sheet', 'Choroidal Thickness');

close all

