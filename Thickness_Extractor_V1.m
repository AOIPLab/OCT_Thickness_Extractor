% Thickness_Extractor_V1
% Created by: Mina Gaffney
% Created on: 2022-09-27

% Description:
%
% Inputs:
% Box/Animal/Follett_VisNeuro_chemical/Data
% LUT with header that contains Image name, and units for scale in row 1,
% col 1 and 2 respectively 
% Col 1 contains scan names, col 2 contains scale
%
% Folder that contains all images 
%
% Folder that contains segmentation spreadsheeets
%
% Optional: Load in folders with segmentations from other observers
%
% Outputs: 

%% Load in the lib folder 
mfilefullpath = matlab.desktop.editor.getActiveFilename; 
mfileshortpath = erase(mfilefullpath, 'Thickness_Extractor_V1.m'); 
lib_path = append(mfileshortpath, 'lib');
addpath(genpath(lib_path));

%% Starting actual script 

% User selects the desired directory:
fullpath = uigetdir;

% Find the first level of folders within the chosen directory:
searchlist_lvl1 = dir(fullpath);

% Fix the DOS era issue with the dir function (loads in the parent
% directories '.' and '..')
searchlist_lvl1 = searchlist_lvl1(~ismember({searchlist_lvl1.name}, {'.', '..'}));

% Only look at folders in first level
searchlist_lvl1 = searchlist_lvl1([searchlist_lvl1.isdir]);

% Ignore "Segmentation_bScans" folder
searchlist_lvl1 = searchlist_lvl1(~ismember({searchlist_lvl1.name}, {'Segmentation_bScans'}));

% Pre- allocate cell array for level 1
level1_folders = cell(length(searchlist_lvl1),1);

% Start level 1 loop
for iii = [1:length(searchlist_lvl1)]
    
    
    % Adding the current folder name to the path: 
    level1_folders{iii} = searchlist_lvl1(iii).name;
    lvl1_paths{iii,1} = [fullpath '\' level1_folders{iii}];

    % Find the second level of folders within the current grader folder:
    searchlist_lvl2 = dir(lvl1_paths{iii,1});
    
    % Fix the DOS era issue with the dir function (loads in the parent
    % directories '.' and '..')
    searchlist_lvl2 = searchlist_lvl2(~ismember({searchlist_lvl2.name}, {'.', '..'}));
    
    
    % EMMA and HANNAH have spreadsheets at different levels... also there's
    % segmentation doccuments in EMMA but not HANNAH
    
    if strcmp(level1_folders{iii,1},'EMMA') == 1 
        % remove the 'segmentation doccuments' folder
        searchlist_lvl2 = searchlist_lvl2(~ismember({searchlist_lvl2.name}, {'Segmentation Documents'}));
        
        % Preallocate cell array for level 2 
        level2_folders = cell(length(searchlist_lvl2),1);
        
        % Start level 2 loop 
        for iv = [1:length(searchlist_lvl2)]
            level2_folders{iv,1} = searchlist_lvl2(iv,1).name;
             
            % Adding the current folder name to the path: 
            %level2_folders{iv} = searchlist_lvl2(iv).name;
            lvl2_paths{iv,1} = [fullpath '\EMMA\' level2_folders{iv} '\Segmentation'];

            % Find the second level of folders within the current grader folder:
            searchlist_lvl3 = dir(lvl2_paths{iv,1});
            
            % Find the spreadsheet within the folder 
            current_spreadsheet = searchlist_lvl3(contains({searchlist_lvl3.name}, {'thickness_data.xlsx'}));
            
            if ~isempty(current_spreadsheet)== 1
            % Load in spreadsheet data
            current_xls_contents = xlsread(fullfile(lvl2_paths{iv,1}, current_spreadsheet.name),'thickness_raw');
            
            % Calculate the min and max for row 1
            current_min = min(current_xls_contents(1,:),[], 2); 
            current_max = max(current_xls_contents(1,:),[], 2);
            elseif isempty(current_spreadsheet)== 1
                current_min = ' ';
                current_max = ' '; 
                
            end
            
            % Saving outputs 
            EMMA_outputs{iv,1} = level2_folders{iv,1};
            EMMA_outputs{iv,2} = current_min;
            EMMA_outputs{iv,3} = current_max;
            
            clear current_spreadsheet current_min current_max
        end
        
    end
        % Make level 2 folders alpha numeric 
    if strcmp(level1_folders{iii,1},'HANNAH') == 1
        % remove the 'segmentation doccuments' folder
        searchlist_lvl2 = searchlist_lvl2(~ismember({searchlist_lvl2.name}, {'Segmentation Documents'}));
        
        % Preallocate cell array for level 2 
        level2_folders = cell(length(searchlist_lvl2),1);
        
        % Start level 2 loop 
        for iv = [1:length(searchlist_lvl2)]
            level2_folders{iv,1} = searchlist_lvl2(iv,1).name;
             
            % Adding the current folder name to the path: 
            %level2_folders{iv} = searchlist_lvl2(iv).name;
            lvl2_paths{iv,1} = [fullpath '\HANNAH\' level2_folders{iv}];

            % Find the second level of folders within the current grader folder:
            searchlist_lvl3 = dir(lvl2_paths{iv,1});
            
            % Find the spreadsheet within the folder 
            current_spreadsheet = searchlist_lvl3(contains({searchlist_lvl3.name}, {'thickness_data.xlsx'}));
            
            if ~isempty(current_spreadsheet)== 1
            % Load in spreadsheet data
            current_xls_contents = xlsread(fullfile(lvl2_paths{iv,1}, current_spreadsheet.name),'thickness_raw');
            
            % Calculate the min and max for row 1
            current_min = min(current_xls_contents(1,:),[], 2); 
            current_max = max(current_xls_contents(1,:),[], 2);
            elseif isempty(current_spreadsheet)== 1
                current_min = ' ';
                current_max = ' '; 
                
            end
            
            % Saving outputs 
            HANNAH_outputs{iv,1} = level2_folders{iv,1};
            HANNAH_outputs{iv,2} = current_min;
            HANNAH_outputs{iv,3} = current_max;
            
            clear current_spreadsheet current_min current_max
        end
    end
     
% Writing outputs to file 
outfile = fullfile(fullpath,['Chemical_OCT_bScans_MATLAB_', datestr(now,'yyyymmdd_HH_MM_SS'), '.xlsx']);
xlswrite(outfile, HANNAH_outputs,'HANNAH');
xlswrite(outfile, EMMA_outputs,'EMMA');


end
