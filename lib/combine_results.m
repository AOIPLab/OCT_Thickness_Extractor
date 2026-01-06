function [results] = combine_results(results, data_dir_mat)
% This function combines data from left and right of the center of the ONH
% data that was determined in earlier steps
   
    for v = 1:length(data_dir_mat)
    
        % combine 0s for total thickness
        results(v).zero_thickness_total = [results(v).zero_vals_left, results(v).zero_vals_right]; % combine the left and right 0 data
        l_r_0_mean = mean(results(v).zero_thickness_total); % get the mean of the left and right 0s thicknesses
        results(v).avg_thickness_val_left(end) = l_r_0_mean; % add mean thickness value to the end of left
        results(v).avg_thickness_val_right(1) = []; % remove first thickness element from right
        results(v).size_of_bin_left(end) = results(v).size_of_bin_left(end) + results(v).size_of_bin_right(1); % add the 0 bin sizes
        results(v).size_of_bin_right(1) = [];% remove first size of bin element from right
        results(v).locations_right(1) = []; % remove first location element from right
        

        % combine left and right values for total (thickness, locations, bin sizes)
        results(v).avg_thickness_total = [results(v).avg_thickness_val_left, results(v).avg_thickness_val_right];
        results(v).locations_total = [results(v).locations_left, results(v).locations_right];
        results(v).size_of_bin_total = [results(v).size_of_bin_left, results(v).size_of_bin_right];
      
    
    end

end