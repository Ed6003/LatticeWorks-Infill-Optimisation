% Set folder path and retrieve CSV files
folderPath = 'D:\TechnicalReport\nTop Homogenisation';
files = getCsvFiles(folderPath);  % Function to get .csv files
fileNames = {files.name};

% Define colour schemes
defineColours();

% Extract numeric values from file names and sort
numericValues = extractNumericValues(fileNames);
sortedFilenames = sortFilenames(fileNames, numericValues);

% Read and store data from CSV files
dataTables = loadDataFromFiles(folderPath, sortedFilenames, numericValues);

% Extract infill densities and stiffness matrices
[infillPercentage, stiffnessData] = extractInfillAndStiffness(dataTables);

% Plot stiffness matrix heatmaps
plotStiffnessHeatmaps(infillPercentage, stiffnessData, sortedFilenames);

% Perform exponential regression on stiffness components
fitResults = performExponentialRegression(infillPercentage, stiffnessData);

% Compare MATLAB results with nTopology results
results = compareResults(stiffnessData);

% Convert results to table and plot
resultsTable = convertResultsToTable(results);
plotOrthotropicProperties(infillPercentage, resultsTable);

% Save figures in specified folder
saveFigures(folderPath, true, 0);