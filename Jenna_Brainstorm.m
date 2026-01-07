% Jenna's Brainstorm
% Created by: Jenna Grieshop
% Date created: 1/7/2025
%
% Attempting to combine Thickness Extractor with Retinal Thickness Analysis
% code to be more elegant

%% Set up environment

clear all
close all
clc

basePath = which('Thickness_Extractor_new.m');

[basePath ] = fileparts(basePath);
path(path,fullfile(basePath,'lib')); % Add our support library to the path.

orientation = 1;

%% User selects data directory, data type determined, data loaded

% user selects data directory
data_dir_path = uigetdir('.','Select directory of data ');


% Get directory of just the .mat files (used for DOCTRAP version)
data_dir = dir(fullfile(data_dir_path, '*.mat'));
if isempty(data_dir)

    % Get directory of just the .csv files (used for generic version)
    data_dir = dir(fullfile(data_dir_path, '*.csv'));
    containsOutput = contains({data_dir.name}, 'LUT');
    data_dir = data_dir(~containsOutput);
    

    % Layers/columns arrainged bottom to top or top to bottom
    list = {'Top to Bottom', 'Bottom to Top'};
    [indx, tf] = listdlg('PromptString', 'Select the orientation of the layer order in the csv sheets.', 'SelectionMode', 'single', 'ListString', list);
    if tf == 1
        orientation = indx; % top to bottom is 1, bottom to top is 2
    else
        % Canceled dialog box - end the program
        return
    end

    % Read files into a variable "OCTdata" 
    for i =1:length(data_dir)
        OCTdata{i} = dlmread( fullfile(data_dir(i).folder, data_dir(i).name) );
    end
    totOCT=length(OCTdata);
    dtype = 2; % for csv general data

else
    % Read files into a variable "OCTdata" (format to match generic entry)
    for i =1:length(data_dir)
        % TODO add warning suppression
        temp = load( fullfile(data_dir(i).folder, data_dir(i).name) );
        temp2 = temp.bScan.CorrectedLayers';
        for j=1:size(temp2,2)
            temp3{j} = [(1:size(temp2,1))', temp2(:,j)];
        end
        OCTdata{i} = [temp3{:}];
    end
    totOCT=length(OCTdata);
    dtype = 1; % for .mat doctrap data
    
end


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

%% Retinal Analysis

for i = 1:totOCT

    % Find the index from the LUT file for the data                     
    for index = 1:length(LUT)
        trigger = find(strcmp(LUT{index,1},(data_dir(i).name)));
        if length(trigger) == 1
            break
        end
    end
    name = data_dir(i).name;

    % Set the axial scale
    vertical_scale = LUT{index, 2};
    
    % set the lateral scale based on the unit
    if unit == (1 || 2) % degrees or microns
        lateral_scale = LUT{index, 3};
    else  % pixels
        lateral_scale = 1;
    end
    
    % load approximate 0 with user selected seed from LUT
    seed_px = LUT{index,4};

    % load the image height from the LUT
    height = LUT{index, 5};

    % extract info from OCT data
    numlayers = size(OCTdata{i},2)/2;

    %% Reading Data Columns 
    % based on orientation
    if orientation == 1
        data = OCTdata{i};
        count = 1;
        for j=1:numlayers
            X_og{j} = data(:, count).* lateral_scale;
            Y_og{j} = data(:, count+1).* vertical_scale;
            count = count+2;
        end 

    else
        data = fliplr(OCTdata{i});
        count = 1;
        for j=1:numlayers
            X_og{j} = data(:, count+1).* lateral_scale;
            Y_og{j} = data(:, count);
            Y_og{j} = height - Y_og{j}; % need to figure out how to get the dimension of the image if generic
            Y_og{j} = Y_og{j}.* vertical_scale;
            count = count+2;
        end
    end

    clear data


    %% Scaling

    % Y_og  = Y_og{i} .* vertical_scale;
    % X_og  = X_og{i} .* lateral_scale;

    % cX=X_og-X_og(seed_px);

    %% Interpolating generic data
    if dtype ==2

        % Defining Column Vectors for interpolation
        % ---------------Define the sampling rate multiplier here------------------
        samplingNum = 6; %This is the amount of interpolation can be changed at will
        % -------------------------------------------------------------------------
        
        % CREATING interpolation vector  for Total thickness ----------------------
        
        % finding the spacing of new sample points? All layers X value min and
        % max the same so just one arbritrarily TODO: check to find actual
        % min and max
        X = ( max(X_og{1}) - min(X_og{1}) ) / ( samplingNum*(length(X_og{1})) );
        N = 0;
        
        % getting new sampling points
        for m = 1 : ( samplingNum*length(X_og{1}) )
            X_TOTAL(m) = X_og{1}(1,1) + N*(X);
            N = N + 1;
        end

        % Interpolating the Layer Thickness
        % Interpolating Y vectors 
        count = 1;
        for j=1:numlayers
            Y{count}  = interp1q(X_og{j}, Y_og{j}, X_TOTAL');
            count = count +1 ;
        end  

        

  
    else
        Y = Y_og;

    end

    cX{j} = X{j} - X{j}(seed_px);

    % This is used for flattening the data
    
    min_Y = min(Y{end});
    [size_bottom,~]   = size(Y{end});
    
    for k = 1:size_bottom
        bottom_shifts(k) = (Y{end}(k) - min_Y);
    end
    
    % shifting all the layers to be relative to the bottom layer
    for j=1:numlayers
        for p = 1:size_bottom 
            Y_shift{j}(p)  = Y{j}(p)  - bottom_shifts(p)';
        end
    end

    %% Finding real 0 from seed point doctrap data (animal only?)

    % look to the left of the seed for the edge (end of 0s in corroidal layer)
    iterationl = 0;
    chor_vall = Y(end-1, seed_px);
    while chor_vall == 0
        iterationl = iterationl + 1;
        chor_vall = Y(end-1, seed_px - iterationl);   
    end
    
    % look to the right of the seed for the edge (end of 0s in corroidal layer)
    iterationr = 0;
    chor_valr = Y(end-1, seed_px);
    while chor_valr == 0
        iterationr = iterationr + 1;
        chor_valr = Y(end-1, seed_px+iterationr);
    end
    
    % find the range of the ONH
    left_range = seed_px - iterationl + 1;
    right_range = seed_px + iterationr - 1;
    
    % find actual 0 based on the range - needs to be an integer so round
    seed_px = round((left_range + right_range) / 2);



    %% Thickness calculations








end




