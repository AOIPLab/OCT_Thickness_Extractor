% Jenna's Brainstorm
% Created by: Jenna Grieshop
% Date created: 1/7/2025
%
% Attempting to combine Thickness Extractor with Retinal Thickness Analysis
% code to be more elegant



%TODO add eye flipping so that they are all same orientation?

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
    

    % % Layers/columns arrainged bottom to top or top to bottom
    % list = {'Top to Bottom', 'Bottom to Top'};
    % [indx, tf] = listdlg('PromptString', 'Select the orientation of the layer order in the csv sheets.', 'SelectionMode', 'single', 'ListString', list);
    % if tf == 1
    %     orientation = indx; % top to bottom is 1, bottom to top is 2
    % else
    %     % Canceled dialog box - end the program
    %     return
    % end

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
list = {'degrees', 'mm', 'pixels'};
[indx_unit,tf1] = listdlg('PromptString','Select Lateral Units', 'SelectionMode', 'single', 'ListString', list);

% set unit text for spacing and window dailog box
if tf1 == 1
    if indx_unit == 1
        unit = 'degrees';
    elseif indx_unit == 2
        unit = 'mm';
    elseif indx_unit == 3
        unit = 'pixels';
    end
end

% set up messages to be displayed to user to set spacing and window in the desired lateral unit
m1 = 'Please enter desired SPACING in';
m2 = 'Please enter the desired WINDOW in';
message1 = [m1, ' ', unit];
message2 = [m2, ' ', unit];

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
    name{i} = data_dir(i).name;

    % Set the axial scale
    vertical_scale = LUT{index, 2};
    
    % set the lateral scale based on the unit
    if indx_unit == 1 % degrees 
        lateral_scale = LUT{index, 3};
        max_unit = 50;
        scaling_factor = (1/lateral_scale);
    elseif indx_unit == 2  % microns
        lateral_scale = LUT{index, 4};
        max_unit = 5;
        scaling_factor = lateral_scale/1000; % in mm rather than microns
    else % pixels
        lateral_scale = 1;
        max_unit = 5000;
        scaling_factor = 1;
    end



    % load approximate 0 with user selected seed from LUT
    seed_px = LUT{index,5};

    % % load the image height from the LUT
    % height = LUT{index, 5};

    % extract info from OCT data
    numlayers = size(OCTdata{i},2)/2;

    %% Reading Data Columns 
    % based on orientation
    if orientation == 1
        data = OCTdata{i};
        count = 1;
        for j=1:numlayers
            if dtype == 1
                X_og{j} = data(:, count).* lateral_scale;
            else
                X_og{j} = data(:, count);
            end
            Y_og{j} = data(:, count+1);
            count = count+2;
        end 

    % else
    %     data = fliplr(OCTdata{i});
    %     count = 1;
    %     for j=1:numlayers
    %         X_og{j} = data(:, count+1).* lateral_scale;
    %         Y_og{j} = data(:, count);
    %         Y_og{j} = height - Y_og{j};
    %         Y_og{j} = Y_og{j};
    %         count = count+2;
    %     end
    end

    clear data


    %% Interpolating generic data
    if dtype ==2

        % Defining Column Vectors for interpolation
        % ---------------Define the sampling rate multiplier here------------------
        % samplingNum = 6; %This is the amount of interpolation can be changed at will
        % -------------------------------------------------------------------------
        
        % CREATING interpolation vector  for Total thickness ----------------------

        % Centering pixels based on seed point
        for q = 1:length(X_og)
            X_centered{q} = X_og{q} - seed_px;
        end
        
        
        % finding the spacing of new sample points? All layers X value min and
        % max the same so just one arbritrarily TODO: check to find actual
        % min and max
        % X = ( max(X_centered{1}) - min(X_centered{1}) ) / ( samplingNum*(length(X_centered{1}))-1 );

        % Assumption/requirement of all the boundaries to start and end
        % with same x coords
        X_TOTAL = min(X_centered{1}):max(X_centered{1});
        
        % % getting new sampling points
        % for m = 1 : ( samplingNum*length(X_centered{1}) )
        %     X_TOTAL(m) = X_centered{1}(1,1) + N*(X);
        %     N = N + 1;
        % end

        % Interpolating the Layer Thickness
        % Interpolating Y vectors 
        count = 1;
        for j=1:numlayers
            Y{count}  = interp1q(X_centered{j}, Y_og{j}, X_TOTAL');
            count = count +1 ;
        end  


        % NOW I think X_TOTAL is the correct pixel eccentricities for the
        % interpolated data - check this - JG

        X_TOTAL_Centered = X_TOTAL * scaling_factor; % now in mm
        X_TOTAL_Centered = X_TOTAL_Centered';
        
    else
        Y = Y_og;
        X_TOTAL = X_og{1};     

    end

    

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


     %% Thickness calculations

    % Total Thickness 
    Ytotal = abs(Y_shift{1} - Y_shift{end});
    % Converting total thickness from pixels to microns
    Ytotal = Ytotal.* vertical_scale;
    Y_layer_thickness{1} = Ytotal;


    c=2;
    for j=1:numlayers-1 
        % thickness between layers
        Y_layer_thickness{c} = abs(Y_shift{j} - Y_shift{j+1});
        Y_layer_thickness{c} = Y_layer_thickness{c}.*vertical_scale;
        c = c+1;
        
    end

    % up to here getting same csv results for thicknesses - no errors - JG
    % up to here getting same .mat results for thicknesses - no errors - JG
    % total thickness appears different as I'm doing last minus 1st
    % segmentation rather than 3-4

    %% Finding real 0 from seed point doctrap data (animal only?) % I think this should be removed for the general version
    if dtype == 1

        % look to the left of the seed for the edge (end of 0s in corroidal layer)
        iterationl = 0;
        chor_vall = Y_layer_thickness{end}(seed_px);
        while chor_vall == 0
            iterationl = iterationl + 1;
            chor_vall = Y_layer_thickness{end}(seed_px - iterationl);   
        end
        
        % look to the right of the seed for the edge (end of 0s in corroidal layer)
        iterationr = 0;
        chor_valr = Y_layer_thickness{end}(seed_px);
        while chor_valr == 0
            iterationr = iterationr + 1;
            chor_valr = Y_layer_thickness{end}(seed_px+iterationr);
        end
        
        % find the range of the ONH
        left_range = seed_px - iterationl + 1;
        right_range = seed_px + iterationr - 1;
        
        % find actual 0 based on the range - needs to be an integer so round
        seed_px = round((left_range + right_range) / 2);
        X_TOTAL_Centered = X_TOTAL - X_TOTAL(seed_px);
    
    else
        % Make make XnewTOTAL work for -x:resolution:+x
        % [dx,numpoints] = size(X_TOTAL); 
        % corrsizex = (round(LUT{index, 9})*LUT{index, 6})/24.46; % assuming 24.46mm eye
        % corrsizexpx = LUT{index, 7};
        % xn = -seed_px+1:0;
        % xp = 1:corrsizexpx-seed_px;
        % comb = [xn, xp];



        % XnewTOTAL = -(corrsizex/2):corrsizex/(numpoints-1):(corrsizex/2); 
        % xl = -(0:corrsizex/(numpoints-1):(corrsizex/2));
        % xr = corrsizex/(numpoints-1):corrsizex/(numpoints-1):(corrsizex/2);
        % xcomb = [fliplr(xl), xr];
        % cb = find(xcomb == 0);


        % X_TOTAL_Centered = XnewTOTAL/1000; %Put the size into 'mm' for later
        % X_TOTAL_Centered = XnewTOTAL';
        
    end

    %% Write All Thickness
    LayersFileOut = strcat([name{i}(1:end-4),'_Retinal_Layers_Thickness.csv']);
    first = 1;
    for j=1:numlayers-1
        if first == 1
            Layersdata = [X_TOTAL_Centered Y_layer_thickness{j}'];
            first = 0;
        else
            Layersdata = [Layersdata Y_layer_thickness{j}'];
        end
    end
    
    %To Write excel file with Sampling vector and thickness values together
    writematrix(Layersdata,LayersFileOut);
    
    clear Ytotal Y_layer_thickness XnewTOTAL Y XnewTOTAL Y_shift

    %% Get bins 

    if dtype == 1

        size_of_matrix_in_unit = size(Layersdata, 1) * scaling_factor;
        sampling_points = -size_of_matrix_in_unit/2:spacing:size_of_matrix_in_unit/2; 
        px_to_unit_n = -(0:scaling_factor:max_unit);
        px_to_unit_p = scaling_factor:scaling_factor:max_unit;
        combined = [fliplr(px_to_unit_n), px_to_unit_p];
        center_bin = find(combined == 0);
    
        right = size(Layersdata, 1) - seed_px;
        eccentricity = combined(center_bin - seed_px+1:center_bin + right);
        Layersdata(:,1) = eccentricity;
    
        % sampling_points = min( Layersdata(:,1)):spacing:max( Layersdata(:,1)+spacing); % Adding spacing to the endpoint ensures that I include all possible points
        sampling_pointsl = round(0:-spacing:min( Layersdata(:,1)-spacing), 4);
        sampling_pointsr = round(0:spacing:max( Layersdata(:,1)+spacing), 4);
        sampling_points = [fliplr(sampling_pointsl) sampling_pointsr(2:end)];
    
        edges = zeros(size(sampling_points,2)+1, 1);
        edges(1,1) = sampling_points(1) - (window/2);
        edges(2:end) = sampling_points + (window/2);

    else
        size_of_matrix_in_unit = size(Layersdata, 1) * scaling_factor;
        % sampling_points = -size_of_matrix_in_unit/2:spacing:size_of_matrix_in_unit/2+spacing; % adding spacing to the endpoint ensures that included all possible points
        % px_to_unit_n = -(0:scaling_factor:max_unit);
        % px_to_unit_p = scaling_factor:scaling_factor:max_unit;
        % combined = [fliplr(px_to_unit_n), px_to_unit_p];
        % center_bin = find(combined == 0);
    
        % right = size(Layersdata, 1) - seed_px;
        % eccentricity = combined(center_bin - seed_px+1:center_bin + right);
        % Layersdata(:,1) = eccentricity;
    
        % sampling_points = min( Layersdata(:,1)):spacing:max( Layersdata(:,1)+spacing); % Adding spacing to the endpoint ensures that I include all possible points
        sampling_pointsl = round(0:-spacing:min( Layersdata(:,1)-spacing),4);
        sampling_pointsr = round(0:spacing:max( Layersdata(:,1)+spacing),4);
        sampling_points = [fliplr(sampling_pointsl) sampling_pointsr(2:end)];
    
        edges = zeros(size(sampling_points,2)+1, 1);
        edges(1,1) = sampling_points(1) - (window/2);
        edges(2:end) = sampling_points + (window/2);
    end
  

    % edges=edges';
    [N, bin] = histc(Layersdata(:,1),edges); 
    Layersdata(:,end+1)=bin;

    filename = [name{i}(1:end-4), '_stats', '.csv'];

    % initialize table
    stats = zeros(length(edges)-1, numlayers*2-1);

    for r = 1:length(edges)-1
        

        [rows] = find((Layersdata(:,end)) == r);

        count = 1;
        first = 1;
        if isempty(rows)

            %insert zero values
            X_value = sampling_points(r);
            for k=2:numlayers
                average{count} = 0;
                % stdev{count} = 0;

                if first == 1
                    stats(r, 1:3) = [X_value average{count} 0];
                    first = 0;
                else
                    stats(r, 2*k-2:2*k-1) = [average{count} 0];
                end
                count = count+1;
            end
        else
            %compute averages
            X_value = sampling_points(r);
            for k=2:numlayers
                average{count} = mean(Layersdata(rows(1):rows(end),k));
                % stdev{count} = std(Layersdata(rows(1):rows(end),k));

                if first == 1
                    stats(r, 1:3) = [X_value average{count} length(Layersdata(rows(1):rows(end),k))];
                    first = 0;
                else
                    stats(r, 2*k-2:2*k-1) = [average{count} length(Layersdata(rows(1):rows(end),k))];
                end
                count = count+1;
            end
            
        end      

    end
    compiled_data{i} = stats;

    % need to average the bins across subjects now
    output_fname = strcat(name{i}(1:end-4), '_','OCT_Thickness_Results_spacing_', spacing_str{1}, '_window_', window_str{1}, '_', unit, '_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.xlsx');        

    
    header_individual = {'Eccentricity'; 'Datapoints Contributed'; 'Average'};
    for a = 1:numlayers-1
        data_compiled_i = [num2cell(stats(:,1)'); num2cell(stats(:, a*2+1))'; num2cell(stats(:, a*2))';];
        data_compiled_individual = [header_individual, data_compiled_i];
        sheet_name = num2str(a);
        writecell(data_compiled_individual, fullfile(LUTpath,output_fname), 'Sheet', a);
    end
        



end

%% Save Summary Data

% initialize min and max location
min_loc = min(compiled_data{1,1}(:,1));
max_loc = max(compiled_data{1,1}(:,1));

for v = 1:length(data_dir)
    if min_loc > min(compiled_data{1,v}(:,1))
        min_loc = min(compiled_data{1,v}(:,1));
    end
    if max_loc < max(compiled_data{1,v}(:,1))
        max_loc = max(compiled_data{1,v}(:,1));
    end
end

% Get the step size betweeen the bin centers
all_bin_locs = round(min_loc:spacing:max_loc,4);

for p = 1:numlayers-1
    results_values{p} = zeros(length(all_bin_locs),length(data_dir));
    results_sum{p} = zeros(length(all_bin_locs),1);
    results_contributed{p} = zeros(length(all_bin_locs),1);
end
%%



% Find which subjects have data at each possible location
for z = 1:length(all_bin_locs')
    for j = 1:length(data_dir)
        max_size_of_bin = max(compiled_data{1,j}(:,3));   
        for k = 1:length(compiled_data{1,j}(:,3))
            if compiled_data{1,j}(k,1) == all_bin_locs(z) % there's an issue with the eccentricities having very very tiny decimals that are making them unequal
                if compiled_data{1,j}(k,3) >= (max_size_of_bin - 1) % make sure to only use full bins
                    for p = 1:numlayers-1

                        results_values{p}(z,j) = compiled_data{1,j}(k,p*2);
                        results_sum{p}(z) = results_sum{p}(z) + compiled_data{1,j}(k,p*2);
                        results_contributed{p}(z) = results_contributed{p}(z) + 1;

                    end
                    
                    break
                end
            end
        end    
    end
end



output_fname_all = strcat('OCT_Thickness_Summary_Results_spacing_', spacing_str{1}, '_window_', window_str{1}, '_', unit, '_', string(datetime('now','TimeZone','local','Format','yyyyMMdd')), '.xlsx');        
header_summary = ['Eccentricity'; 'Datapoints Contributed'; 'Average'; name(:)];
for a = 1:numlayers-1
    data_compiled_s = [num2cell(all_bin_locs); num2cell(results_contributed{a})'; num2cell((results_sum{a}./results_contributed{a}))'; num2cell(results_values{a})'];
    data_compiled_summary = [header_summary, data_compiled_s];
    sheet_name = num2str(a);
    writecell(data_compiled_summary, fullfile(LUTpath,output_fname_all), 'Sheet', a);
end
