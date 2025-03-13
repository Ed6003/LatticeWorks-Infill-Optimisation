% Clear environment and close all figures
clear all; close all; clc;

% Set default folder for results
defaultFolder = 'D:\TechnicalReport\LatticeProperties_LinearFinal';
matPath = fullfile(defaultFolder, '!Summary', 'simulation_results_summary.mat');

% Check if summary results file exists
if fileExists(matPath)
    load(matPath);
else
    % Initialise variables for processing
    count = 1;
    linspace_steps = 64; 
    nPoints = 64; 
    original_points = linspace(1.3, -1.3, linspace_steps);
    selected_idx = round(linspace(1, linspace_steps, nPoints));
    exclusions = {'f', 'v', 'c', 'F', 'V', 'meshOutput', 'abaqusData', 'E_effectiveStrain', 'E_effectiveStress'};

    % Loop through selected points
    for i = original_points(selected_idx)
        savePath = fullfile(defaultFolder, sprintf('%.5g_Lattice_Density', i));
        load(matlabPath);

        % Get relevant fields, excluding specified ones
        expectedFields = setdiff(fieldnames(resultStruct), exclusions);

        % Preallocate summary structure on first iteration
        if count == 1
            template = cell2struct(cell(size(expectedFields)), expectedFields, 1);
            summaryStruct = repmat(template, 1, linspace_steps);
        end

        % Populate summary structure with existing fields
        for field in expectedFields
            summaryStruct(count).(field) = getFieldValue(resultStruct, field);
        end
        
        disp(count); 
        count++;
    end

    % Save summary results
    save(matPath, 'summaryStruct', '-v7.3'); % Use version for large sizes
end

% Define colours for plots
blue = [0, 0.4470, 0.7410];
cyan = [0.3010, 0.7450, 0.9330];

% ---- Plotting Section ----

% Figure 1: Poisson's Ratio
plotPoissonsRatio(summaryStruct, blue, cyan);

% Figure 2: Young's Modulus Fit
fitYoungsModulus(summaryStruct, blue, cyan);

% Figure 3: Young's Moduli R^2 Fit
plotYoungsModuliR2(summaryStruct, blue, cyan);

% Figure 4: Distorted Elements
plotDistortedElements(summaryStruct, blue);

% Figure 5: Infill Percentage vs Level-Set
plotLevelSet(summaryStruct, blue);

% Figure 6: Simulation Time
plotSimulationTime(summaryStruct, blue);

% Figure 7: CPU Load Moving Average
plotCPULoadMovingAverage();

% Figure 8: Number of Tracked Nodes
plotTrackedNodes(summaryStruct, blue, cyan);

% Figure 9: Stress at Each Infill Percentage
plotStress(summaryStruct, blue);

% Figure 10: Poisson's Ratio with Uncertainty
plotPoissonsRatioWithUncertainty(summaryStruct, blue, cyan);

% Figure 11: Young's Moduli with Uncertainty
plotYoungsModuliWithUncertainty(summaryStruct, blue, cyan);

% Save figures in summary folder
saveFiguresInSummaryFolder(defaultFolder);