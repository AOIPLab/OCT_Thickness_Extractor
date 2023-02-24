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

%% Loading in LUT file, OCT image folder, Observer segmentation spreasheets 

% User selects the LUT:
[LUTfile, LUTpath] = uigetfile('.xlsx','Select the LUT');
LUT = readcell(fullfile(LUTpath,LUTfile));

% User selects the desired directory containing all images:
img_fullpath = uigetdir('.','Select directory containing images');

% Find the files within the chosen directory:
img_folder = dir(img_fullpath);

% Fix the DOS era issue with the dir function (loads in the parent
% directories '.' and '..')
img_folder = img_folder(~ismember({img_folder.name}, {'.', '..'}));

% Pre- allocate cell array for all images
OCT_images = cell(length(img_folder),1);

% Select segmentations from observer 1
obs1_fullpath = uigetdir('.','Select directory containing observer 1 segmentation spreadsheets');

% Find the files within the chosen directory:
obs1_folder = dir(obs1_fullpath);

% Fix the DOS era issue with the dir function (loads in the parent
% directories '.' and '..')
obs1_folder = obs1_folder(~ismember({obs1_folder.name}, {'.', '..'}));

% Pre- allocate cell array for all segmentations
obs1_segmentation_folders = cell(length(obs1_folder),2);
obs1_segmentation_folders = struct2cell(obs1_folder);


% Prompt for pop-up asking if you have more than one observer
observers = questdlg('Add another observer', ...
      'Add observer', ...
      'YES', 'NO', 'NO');
  
  if strcmpi(observers, 'NO')
  elseif strcmpi(observers, 'YES')
    % Select segmentations from observer 1
    obs2_fullpath = uigetdir('.','Select directory containing observer 2 segmentation spreadsheets');

    % Find the files within the chosen directory:
    obs2_folder = dir(obs2_fullpath);
    
    % Fix the DOS era issue with the dir function (loads in the parent
    % directories '.' and '..')
    obs2_folder = obs2_folder(~ismember({obs2_folder.name}, {'.', '..'}));

    % Pre- allocate cell array for all images
    obs2_segmentation_folders = cell(length(obs2_folder),1);
    
    % 
  end

% Start level 1 loop
for iii = [1:length(img_folder)]
    
    
    % Adding the current folder name to the path: 
    level1_folders{iii} = img_folder(iii).name;
    lvl1_paths{iii,1} = [img_fullpath '\' level1_folders{iii}];

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
            lvl2_paths{iv,1} = [img_fullpath '\EMMA\' level2_folders{iv} '\Segmentation'];

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
            lvl2_paths{iv,1} = [img_fullpath '\HANNAH\' level2_folders{iv}];

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
outfile = fullfile(img_fullpath,['Chemical_OCT_bScans_MATLAB_', datestr(now,'yyyymmdd_HH_MM_SS'), '.xlsx']);
%xlswrite(outfile, HANNAH_outputs,'HANNAH');
%xlswrite(outfile, EMMA_outputs,'EMMA');


end
