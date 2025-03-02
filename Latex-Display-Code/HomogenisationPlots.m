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
tables = cell(length(sortedFilenames), 1);
for i = 1:length(sortedFilenames)
    filePath = fullfile(folderPath, sortedFilenames{i});
    tables{i} = readtable(filePath);
    tables{i,2} = numericValues(sortedIndices(i)) * 100; % convert to percentage
end

densities = vertcat(tables{:,2});

%(*@\codesubsection{Stiffness Matrix Heatmaps}{homogenisation-stiffness-matrix-heatmaps}@*)
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
    title(sprintf('%.1f%% Density', densities(idx)));
    set(gca, 'XTick', 1:6, 'YTick', 1:6);
end
sgtitle('Stiffness Matrices (C)');

%(*@\codesubsection{Exponential Regression}{homogenisation-stiffness-exponential-regression}@*)
cFigure;
hold on;

% preallocate
fitResults = cell(6,1); % stores (a, b, R^2)

for diagIdx = 1:6
    x = densities(:);
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

% compile results into table
componentNames = cell(6,1); aValues = zeros(6,1);
bValues = zeros(6,1); r2Values = zeros(6,1);

for diagIdx = 1:6
    componentNames{diagIdx} = sprintf('C_{%d%d}', diagIdx, diagIdx);
    if isstruct(fitResults{diagIdx})
        res = fitResults{diagIdx};
        aValues(diagIdx) = res.a;
        bValues(diagIdx) = res.b;
        r2Values(diagIdx) = res.Rsquare;
    else
        aValues(diagIdx) = NaN;
        bValues(diagIdx) = NaN;
        r2Values(diagIdx) = NaN;
    end
end

resultsTable = table(componentNames, aValues, bValues, r2Values, 'VariableNames', {'Component', 'a', 'b', 'R2'});
disp(resultsTable);

%(*@\codesubsection{Save Figures}{homogenisation-save-figures}@*)
saveFigures(folderPath,true,0);