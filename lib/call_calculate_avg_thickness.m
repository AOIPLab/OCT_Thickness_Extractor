%% Call_calculate_avg_thickness
function [results, results_tot, results_in, results_out, results_chor] = call_calculate_avg_thickness(dtype, folder, spacing, window, LUT, unit)
    
    
    % initializing results structs
    results = struct([]);
    results_tot = struct([]);
    results_in = struct([]);
    results_out = struct([]);
    results_chor = struct([]);
    
    %% get path, load segmentation, find lateral scale, set user seed
    
    % this for loop goes through each image folder in the grader folder
    for i = (1:length(folder))
        
        % get the correct path of where the image and segmentation is
        folder_path = [folder(i).folder];
        image_path = [folder(i).folder '\' folder(i).name(1:end-4) '.tif'];
        segmentation_path = [folder(i).folder '\' folder(i).name];
        
        %load the segmentation
        seg = load(segmentation_path);  
        segmentation = seg.bScan.CorrectedLayers;
        
        % add file names to results structs
        % results_tot(i).file_name = folder(i).name;
        % results_in(i).file_name = folder(i).name;
        % results_out(i).file_name = folder(i).name;
        % results_chor(i).file_name = folder(i).name;

        % find the index for the specific scan from
        % the LUT file
        for index = 1:length(LUT)
            trigger = find(strcmp(LUT{index,1},(folder(i).name)));
            if length(trigger) == 1
                break
            end
        end

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

        % call calc thickness function
       
        % this is where we'd have to know the flag of what kind of data we
        % are working with - JG generalize
        if dtype == 1
            scaled_segmentation = funCalcThickness(segmentation, vertical_scale, lateral_scale, seed_px);
        else
            % instead of calling funCalcThickness need to scale the data
            % and then find the thicknesses for each layer...maybe could be
            % done in the same way though and could use the
            % funCalcThickness
        end
        
        %%  find actual 0 from seed reference point, split up the segmetation based on left and right of 0, call calculation function
        
        % look to the left of the seed for the edge (end of 0s in corroidal layer)
        iterationl = 0;
        chor_vall = scaled_segmentation(5, seed_px);
        while chor_vall == 0
            iterationl = iterationl + 1;
            chor_vall = scaled_segmentation(5, seed_px - iterationl);   
        end
        
        % look to the right of the seed for the edge (end of 0s in corroidal layer)
        iterationr = 0;
        chor_valr = scaled_segmentation(5, seed_px);
        while chor_valr == 0
            iterationr = iterationr + 1;
            chor_valr = scaled_segmentation(5, seed_px+iterationr);
        end
        
        % find the range of the ONH
        left_range = seed_px - iterationl + 1;
        right_range = seed_px + iterationr - 1;
        
        % find actual 0 based on the range - needs to be an integer so round
        seed_px = round((left_range + right_range) / 2);

        % the segmentation matrix split up into left and right of the seed
        left_of_seed = scaled_segmentation(:, 1:seed_px);
        left_of_seed = flip(left_of_seed, 2); % flip the matrix so that the left of the seed can be read left to right
        right_of_seed = scaled_segmentation(:, seed_px+1:end);
        
        
        % call 2nd function for the left and the right side of 0
        [results, results_tot, results_in, results_out, results_chor] = calculate_avg_thickness(window, spacing, left_of_seed, 'left', folder, i, results, results_tot, results_in, results_out, results_chor, lateral_scale);
        [results, results_tot, results_in, results_out, results_chor] = calculate_avg_thickness(window, spacing, right_of_seed, 'right', folder, i, results, results_tot, results_in, results_out, results_chor, lateral_scale);

    end
end