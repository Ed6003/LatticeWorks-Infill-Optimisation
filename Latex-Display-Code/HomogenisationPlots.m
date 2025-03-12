folderPath = 'D:\TechnicalReport\nTop Homogenisation'; %(*@\draftcomment{the prefix for every sub-heading in this script is "homogenisation"}@*)
files = dir(fullfile(folderPath, '*.csv'));
fileNames = {files.name};

blue   = [0, 0.4470, 0.7410]; % "#0072BD"
cyan  = [0.3010, 0.7450, 0.9330]; % "#4DBEEE"
red = [0.8500, 0.3250, 0.0980]; % "#D95319"
yellow = [0.9290, 0.6940, 0.1250]; % "#EDB120"

% extract numeric values from fileNames
numericValues = zeros(length(fileNames), 1);
for i = 1:length(fileNames)
    [~, name, ~] = fileparts(fileNames{i});  % remove .csv extension
    numericValues(i) = str2double(name);
end

% sort based on number
[~, sortedIndices] = sort(numericValues);
sortedFilenames = fileNames(sortedIndices);

% data from CSV files in sorted order
dataTables = cell(length(sortedFilenames), 1);
for i = 1:length(sortedFilenames)
    filePath = fullfile(folderPath, sortedFilenames{i});
    dataTables{i} = readtable(filePath);
    dataTables{i,2} = numericValues(sortedIndices(i)) * 100; % convert to percentage
end

%(*@\codesubsection{Extract Infill Densities}{homogenisation-extract-infill-densities}@*)
% Extract infill densities and stiffness matrices
numFiles = length(dataTables);
infillPercentage = zeros(numFiles, 1);
stiffnessData = zeros(numFiles, 6, 6); % 3D array: [density, row, col]

for i = 1:numFiles
    infillPercentage(i) = dataTables{i,2};        % Get density from second column
    currentTable = dataTables{i,1};       % Get 6x6 table from first column
    stiffnessData(i, :, :) = currentTable{:,:}; % Convert table to numeric array
end

%(*@\codesubsection{Stiffness Matrix Heatmaps}{homogenisation-stiffness-matrix-heatmaps}@*)
numFiles = numel(sortedFilenames);
numSteps = 6;
selected_indices = round(linspace(1, numFiles, min(numFiles,numSteps)));

nCols = ceil(sqrt(numSteps));
nRows = ceil(numSteps / nCols);

cFigure;
for j = 1:length(selected_indices)
    idx = selected_indices(j);
    subplot(nRows,nCols,j);

    imagesc(squeeze(stiffnessData(idx, :, :))); % create heatmap
    colorbar;
    axis square;
    xlabel('Column Index');
    ylabel('Row Index');
    title(sprintf('%.1f%% Density', infillPercentage(idx)));
    set(gca, 'XTick', 1:6, 'YTick', 1:6);
end
sgtitle('Stiffness Matrices (C)');

%(*@\codesubsection{Exponential Regression}{homogenisation-stiffness-exponential-regression}@*)
cFigure;
hold on;

% preallocate
fitResults = cell(6,1); % stores (a, b, R^2)

for diagIdx = 1:6
    x = infillPercentage(:);
    y = squeeze(stiffnessData(:, diagIdx, diagIdx));
    y = y(:);

    % original data
    if ismember(diagIdx, [1, 2, 3])
        colour = blue;
    else
        colour = cyan;
    end

    plot(x, y, 'LineWidth', 1.5, 'DisplayName', sprintf('C_{%d%d}', diagIdx, diagIdx),'Color',colour);


    % using custom exponential model y = a*x*exp(b*x) to ensure 0 crossing
    ft = fittype('a*x*exp(b*x)', 'independent', 'x', 'dependent', 'y');

    % linear initial guess (estimating to ensure convergence)
    valid = (x > 0) & (y > 0);
    if sum(valid) >= 2
        xValid = x(valid);
        yValid = y(valid);
        z = log(yValid) - log(xValid);
        X = [ones(size(xValid)), xValid];
        coeffs = X \ z;
        logAInit = coeffs(1);
        bInit = coeffs(2);
        aInit = exp(logAInit);
    end

    [fitResult, gof] = fit(x, y, ft, 'StartPoint', [aInit, bInit]);

    fitResults{diagIdx} = struct('a', fitResult.a, 'b', fitResult.b, 'Rsquare', gof.rsquare);

    % plot fitted curves
    if ismember(diagIdx, [1, 2, 3])
        colour = red;
    else
        colour = yellow;
    end

    plot(x, fitResult(x), '--', 'LineWidth', 1.5, 'DisplayName', sprintf('C_{%d%d} Fit', diagIdx, diagIdx),'Color',colour);
end

xlabel('Infill Density (%)');
ylabel('Stiffness Component Value');
title('Diagonal Stiffness Components against Infill Density');
legend('Location', 'best');
grid on;

%(*@\codesubsection{Matlab and nTopology result comparison}{homogenisation-matlab-and-ntopology-result-comparison}@*)
numFiles = length(stiffnessData);
results = zeros(numFiles, 9);

for i = 1:numFiles
    C = squeeze(stiffnessData(i,:,:)); % Extract 6x6 from 1x6x6
    S = inv(C); 

    % Extract engineering constants
    E1 = 1 / S(1,1);
    E2 = 1 / S(2,2);
    E3 = 1 / S(3,3);
    
    nu12 = -S(2,1) / S(1,1);
    nu13 = -S(3,1) / S(1,1);
    nu23 = -S(3,2) / S(2,2);
    
    G12 = 1 / S(4,4);
    G23 = 1 / S(5,5);
    G31 = 1 / S(6,6);

    % Store results
    results(i, :) = [E1, E2, E3, nu12, nu13, nu23, G12, G23, G31];
end

% Convert to table and display
propertyNames = {'E_x', 'E_y', 'E_z', 'v_{xy}', 'v_{xz}', 'v_{yz}', 'G_{xy}', 'G_{yz}', 'G_{xz}'};
resultsTable = array2table(results, 'VariableNames', propertyNames);

[numFiles, numProperties] = size(resultsTable);

% convert table to array for plotting
data = table2array(resultsTable);
propertyNames = resultsTable.Properties.VariableNames;

% subplots
figure;
hold on;
for i = 1:numProperties
    subplot(3,3,i); % Arrange in a 3x3 grid
    hold on;
    plot(infillPercentage, data(:,i), '-', 'LineWidth', 1.5, 'MarkerSize', 8);
    switch i
        case 1
            plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xy_mean).*1e6,'Color',cyan, 'LineWidth', 1.5)
        case 2
            plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xz_mean).*1e6,'Color',cyan, 'LineWidth', 1.5)
        case 3
            plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xz_mean).*1e6,'Color',cyan, 'LineWidth', 1.5)
        case 4
            plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xy_mean),'Color',cyan, 'LineWidth', 1.5)
        case 5
            plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xz_mean),'Color',cyan, 'LineWidth', 1.5)
        case 6
            plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xy_mean),'Color',cyan, 'LineWidth', 1.5)
    end
    title(propertyNames{i}, 'FontSize', 9);
    xlabel('Infill Percentage (%)');
    ylabel(propertyNames{i});
    grid on;
end
sgtitle('Orthotropic Properties from 6x6 Stiffness Matrix, nTop and Matlab Method Comparison');

%(*@\codesubsection{Save Figures}{homogenisation-save-figures}@*)
saveFigures(folderPath,true,0);