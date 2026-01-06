%% Calculate - Does thickness calculations
function [results, results_tot, results_in, results_out, results_cor] = calculate_avg_thickness(window, spacing, matrix, LorR, folder, i, results, results_tot, results_in, results_out, results_cor, lateral_scale)
 
    % at this point nothing is in the individual structs - JG

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
    first_storejg = 1;
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
        
        for n=2:5 % NEED TO DO THIS BY HOW MANY LAYERS NOT HARD CODED - JG
            values(n-1, window_count) = mean(matrix(n, (index_min:index_max)));
            if first_storejg
                values_zero{n-1} = matrix(n, (index_min:index_max));
            end
        end
        first_storejg = 0;

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
            % add for loop here to go through for all layers - JG
            for n=2:5 % NEED TO DO THIS BY HOW MANY LAYERS NOT HARD CODED - JG
                results(i,n-1).file_name = folder(i).name;
                results(i,n-1).zero_vals_right = values_zero(n-1);
                results(i,n-1).avg_thickness_val_right = values(n-1);
                results(i,n-1).locations_right = sampling_points;
                results(i,n-1).size_of_bin_right = bin_size;
            end

            % store thicknesses
            results_tot(i).zero_vals_right = values_total_zero;
            results_tot(i).avg_thickness_val_right = values_total;

            results_in(i).zero_vals_right = values_inner_zero;
            results_in(i).avg_thickness_val_right = values_inner;

            results_out(i).zero_vals_right = values_outer_zero;
            results_out(i).avg_thickness_val_right = values_outer;

            results_cor(i).zero_vals_right = values_choroidal_zero;
            results_cor(i).avg_thickness_val_right = values_choroidal;
    
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

            % add for loop here to go through for all layers - JG
            for n=2:5 % NEED TO DO THIS BY HOW MANY LAYERS NOT HARD CODED - JG
                results(i,n-1).file_name = folder(i).name;
                results(i,n-1).zero_vals_left = flip(values_zero(n-1));
                results(i,n-1).avg_thickness_val_left = flip(values(n-1));
                results(i,n-1).locations_left = -flip(sampling_points);
                results(i,n-1).size_of_bin_left = flip(bin_size);
            end

            % store thicknesses - flip left results back to read left to right
            results_tot(i).zero_vals_left = flip(values_total_zero);
            results_tot(i).avg_thickness_val_left = flip(values_total,2);
            results_in(i).zero_vals_left = flip(values_inner_zero);
            results_in(i).avg_thickness_val_left = flip(values_inner,2);
            results_out(i).zero_vals_left = flip(values_outer_zero);
            results_out(i).avg_thickness_val_left = flip(values_outer,2);
            results_cor(i).zero_vals_left = flip(values_choroidal_zero);
            results_cor(i).avg_thickness_val_left = flip(values_choroidal,2);
    
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
