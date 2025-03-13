cd('C:\Users\rusco\OneDrive - University of Warwick\Admin\Archives\Documents\GitHub\LatticeWorks');
checkSubfolders = true;

saveFolder = 'C:\Users\rusco\OneDrive - University of Warwick\Admin\Archives\Documents\GitHub\LatticeWorks\LogsActual';
masterFile = fullfile(saveFolder, 'master.csv');

logsDir = 'C:\Users\rusco\OneDrive - University of Warwick\Admin\Archives\Documents\GitHub\LatticeWorks\LogsActual';

% ensure the Logs directory exists
if ~isfolder(logsDir)
    error('The Logs directory does not exist.');
end

% get a list of all subfolders in the Logs directory
subfolders = dir(logsDir);
subfolders = subfolders([subfolders.isdir]); % keep only directories
subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % exclude . and ..

if isempty(subfolders)
    checkSubfolders = false;
end

% find the most recently modified subfolder
[~, sortIdx] = sort([subfolders.datenum]); % Sort by creation time
subfolders = subfolders(sortIdx);

if checkSubfolders % creating the csv files from the blg files
    for i = 1:length(subfolders)
        recentSubfolder = subfolders(i).name;
        recentSubfolderPath = fullfile(logsDir, recentSubfolder);

        % locate the .blg file in the most recent subfolder
        blgFiles = dir(fullfile(recentSubfolderPath, '*.blg'));

        % check if a .blg file exists
        if isempty(blgFiles)
            error('No .blg file found in the most recent subfolder, check folder permissions. (Admin only by default)');
        end

        % assume there's only one .blg file; if multiple, process the first one
        blgFile = fullfile(recentSubfolderPath, blgFiles(1).name);

        % define the output .csv file path
        csvName = sprintf('%.0f', posixtime(datetime('now', 'TimeZone', 'UTC')));
        csvFile = fullfile(logsDir, [csvName, '.csv']);

        % convert .blg to .csv using the relog command
        relogCmd = sprintf('relog "%s" -f csv -o "%s"', blgFile, csvFile);
        system(relogCmd);

        pause(1);

        try
            rmdir(recentSubfolderPath, 's');
        catch
            warning('Failed to delete the folder: %s\nError: %s', recentSubfolderPath, ME.message);
        end
    end
end


% get list of CSV files sorted by creation date
fileList = dir(fullfile(saveFolder, '*.csv'));
fileList = fileList(~strcmp({fileList.name}, 'master.csv'));
[~, sortIdx] = sort([fileList.datenum]); % Sort by creation time
fileList = fileList(sortIdx);

% loop through each file in order of creation and append to master.csv
for i = 1:length(fileList)
    filePath = fullfile(saveFolder, fileList(i).name);

    % Read the CSV file
    data = readtable(filePath);

    colIdx = contains(data.Properties.VariableNames, 'Total') & contains(data.Properties.VariableNames, 'ProcessorTime') | contains(data.Properties.VariableNames, 'GMTStandardTime');

    columnName = data.Properties.VariableNames(colIdx);
    data = data(:, columnName);

    data.Properties.VariableNames(1) = "Time";
    data.Properties.VariableNames(2) = "Processor Time %";

    numericVars = varfun(@isnumeric, data, 'OutputFormat', 'uniform');
    data(:, numericVars) = fillmissing(data(:, numericVars), 'constant', 0);

    % append to master.csv (create if it doesn't exist)
    if exist(masterFile, 'file')
        writetable(data, masterFile, 'WriteMode', 'append');
    else
        writetable(data, masterFile);
    end

    delete(filePath);
end

disp('Processing complete. All files merged and deleted.');