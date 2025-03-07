clear; close all; clc; %(*@\draftcomment{the prefix for every sub-heading in this script is "summarystruct"}@*) %(*@\codesubsection{Merge Results}{summarystruct-merge-structs}@*)
defaultFolder = fullfile('D:','TechnicalReport','LatticeProperties_LinearFinal');
matPath = fullfile(defaultFolder,'!Summary','simulation_results_summary.mat');

if exist(matPath,"file") ~= 0
    load(matPath);
else
    count = 1;

    linspace_steps = 64;
    nPoints = 64; % resample if required, if nto set to = linspace_steps
    original_points = linspace(1.3, -1.3, linspace_steps);
    selected_idx = round(linspace(1, linspace_steps, nPoints));

    exclusions = {'f', 'v', 'c', 'F', 'V', 'meshOutput', 'abaqusData', 'E_effectiveStrain', 'E_effectiveStress'};

    for i = original_points(selected_idx)

        defaultFolder = fullfile('D:','TechnicalReport','LatticeProperties'); % changed to external hard drive
        savePath=fullfile(defaultFolder,sprintf('%.5g_Lattice_Density',i));
        matlabPath=fullfile(savePath,'simulation_results.mat');

        load(matlabPath);

        expectedFields = fieldnames(resultStruct);
        expectedFields = expectedFields(~ismember(expectedFields, exclusions));

        if count == 1
            % pre-allocate
            template = cell2struct(cell(size(expectedFields)), expectedFields, 1);
            summaryStruct = repmat(template, 1, linspace_steps);
        end

        % check if it exists, if not leave empty
        for j = 1:numel(expectedFields)
            fieldName = expectedFields{j};
            if isfield(resultStruct, fieldName)
                summaryStruct(count).(fieldName) = resultStruct.(fieldName);
            else
                summaryStruct(count).(fieldName) = [];  % Or assign a default value
            end
        end
        count = count + 1;
    end

    matlabPath=fullfile(defaultFolder,'!Summary','simulation_results_summary.mat');
    save(matlabPath,'summaryStruct', '-v7.3') % need to specify .mat version to support larger size above 0.5 GB
end

blue   = [0, 0.4470, 0.7410];  % "#0072BD"
cyan  = [0.3010, 0.7450, 0.9330];  % "#4DBEEE"

%(*@\codesubsection{Figure 1}{summarystruct-figure-1}@*)
cFigure;
hold on; grid on;

% compute error as difference between mean and median
error_xy = abs(vertcat(summaryStruct.poisson_xy_mean) - vertcat(summaryStruct.poisson_xy_median));
error_xz = abs(vertcat(summaryStruct.poisson_xz_mean) - vertcat(summaryStruct.poisson_xz_median));

h1 = plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xy_mean), 'Color', blue, 'LineWidth', 1.5);
h2 = plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xy_median), '--', 'Color', blue, 'LineWidth', 1.5);
h3 = plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xz_mean), 'Color', cyan, 'LineWidth', 1.5);
h4 = plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xz_median), '--', 'Color', cyan, 'LineWidth', 1.5);

% error bars, assuming between mean and median
errorbar(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xy_mean), error_xy, 'o', 'Color', blue, 'CapSize', 5);
errorbar(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xz_mean), error_xz, 'o', 'Color', cyan, 'CapSize', 5);

title("Calculated Poisson's Ratio", 'FontSize', 14);
legend([h1, h2, h3, h4], 'XY Mean', 'XY Median', 'XZ Mean', 'XZ Median', 'Location', 'best');
xlabel("Infill Percentage (%)");
ylabel("Poisson's Ratio");

%(*@\codesubsection{Figure 2}{summarystruct-figure-2}@*)
youngs_modulus = vertcat(summaryStruct.youngs_modulus);

cFigure;
hold on; grid on;
plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xy_mean),'Color',blue, 'LineWidth', 1.5)
plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xy_median), '--','Color',blue, 'LineWidth', 1.5)
plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xz_mean),'Color',cyan, 'LineWidth', 1.5)
plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xz_median), '--','Color',cyan, 'LineWidth', 1.5)

title("Calculated Young's Moduli",'FontSize',14)
legend('XZ Mean', 'XY Mean', 'XZ Median', 'XY Median','Location','best');
xlabel("Infill Percentage (%)")
ylabel("Young's Modulus (MPa)")

%(*@\codesubsection{Figure 3}{summarystruct-figure-3}@*)
cFigure;
hold on; grid on;
plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xy_mean_R_squared),'Color',blue, 'LineWidth', 1.5)
plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xy_median_R_squared), '--','Color',blue, 'LineWidth', 1.5)
plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xz_mean_R_squared),'Color',cyan, 'LineWidth', 1.5)
plot(vertcat(summaryStruct.infill_percentage), vertcat(youngs_modulus.youngs_modulus_xz_median_R_squared), '--','Color',cyan, 'LineWidth', 1.5)

title("Calculated Young's Moduli R^2 Fit",'FontSize',14)
legend('XZ Mean', 'XY Mean', 'XZ Median', 'XY Median','Location','best');
xlabel("Infill Percentage (%)")
ylabel("Young's Modulus R^2 Linear Fit")

%(*@\codesubsection{Figure 4}{summarystruct-figure-4}@*)
cFigure;
hold on; grid on;
plot(vertcat(summaryStruct.infill_percentage), (vertcat(summaryStruct.distortedElements)./vertcat(summaryStruct.numElements)*100),'Color',blue, 'LineWidth', 1.5)

title("Distorted Elements",'FontSize',14)
xlabel("Infill Percentage (%)")
ylabel("Distorted Elements (%)")

%(*@\codesubsection{Figure 5}{summarystruct-figure-5}@*)
cFigure;
hold on; grid on;
plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.level_set),'Color',blue, 'LineWidth', 1.5)

title("Infill Percentage against Level-Set",'FontSize',14)
xlabel("Infill Percentage (%)")
ylabel("Gyroid TPMS Level-Set")

%(*@\codesubsection{Figure 6}{summarystruct-figure-6}@*)
timestampStrings = {summaryStruct.timestamp};
timestamps = datetime(timestampStrings, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

timeDifferences = [NaN, minutes(diff(timestamps))]; % add NaN for the first entry
infillPercentages = vertcat(summaryStruct.infill_percentage);

cFigure;
hold on; grid on;

timeDifferences(timeDifferences > 10) = NaN; % Replace with NaN to use fillmissing
timeDifferencesFilled = fillmissing(timeDifferences, 'linear');

plot(infillPercentages, timeDifferencesFilled,'Color',blue, 'LineWidth', 1.5);

xlabel("Infill Percentage (%)");
ylabel("Time Taken (minutes)");
title("Simulation Time at each Infill Percentage",'FontSize',14);

%(*@\codesubsection{Figure 7}{summarystruct-figure-7}@*)
master_csv = 'C:\Users\rusco\OneDrive - University of Warwick\Admin\Archives\Documents\GitHub\LatticeWorks\LogsActual\master.csv';

data = readtable(master_csv);
data.Time = datetime(data.Time, 'InputFormat', 'dd/MM/yyyy HH:mm:ss.SSS');

cFigure;
hold on; grid on;

time_indices = (1:numel(data.ProcessorTime_)) * 5; % 5 comes from the logging frequency (1 every 5 seconds)
time_indices = time_indices / 60^2; % convert to hours

windowSize = 3600/5; % hourly moving average (accounts for 5s sampling)
smoothed_data = movmean(data.ProcessorTime_, windowSize);
plot(time_indices, smoothed_data,'Color',blue, 'LineWidth', 1)

title("CPU Load Hourly Moving Average",'FontSize',14)
xlabel("Computing Hours")
ylabel("CPU Load (%)")

%(*@\codesubsection{Figure 8}{summarystruct-figure-8}@*)
numEntries = size(summaryStruct,2);
TrackZ = zeros(1,numEntries);
TrackY = zeros(1,numEntries);

for i = 1:numEntries
    TrackZ(i) = numel(summaryStruct(i).bcSets.bcTrackZ);
    TrackY(i) = numel(summaryStruct(i).bcSets.bcTrackY);
end

cFigure;
hold on; grid on;

plot(infillPercentages, TrackY,'Color',blue, 'LineWidth', 1.5);
plot(infillPercentages, TrackZ,'Color',cyan, 'LineWidth', 1.5);

xlabel("Infill Percentage (%)");
ylabel("Tracked Nodes");
title("Tracked Nodes at each Infill Percentage",'FontSize',14);
legend('Y Tracked','Z Tracked','Location','best');

%(*@\codesubsection{Figure 9}{summarystruct-figure-9}@*)
stressSteps = zeros(1,numEntries);

cFigure;
hold on; grid on;
for i = 1:numEntries
    stressSteps(i) = summaryStruct(i).stress(end);
end

plot(infillPercentages, stressSteps,'Color',blue, 'LineWidth', 1.5);

xlabel("Infill Percentage (%)");
ylabel("Stress (MPa)");
title("Stress at each Infill Percentage",'FontSize',14);
% legend('Y Tracked','Z Tracked','Location','best');

%(*@\codesubsection{Figure 10}{summarystruct-figure-10}@*)
% plot error avg
figure; 
hold on; grid on;
plot(vertcat(summaryStruct.infill_percentage), error_xy, 'Color', blue, 'LineWidth', 1.5);
plot(vertcat(summaryStruct.infill_percentage), error_xz, 'Color', cyan, 'LineWidth', 1.5);

title("Average Error in Poisson's Ratio", 'FontSize', 14);
xlabel("Infill Percentage (%)");
ylabel("Average Absolute Error");
legend({'XY Error', 'XZ Error'}, 'Location', 'best');

savePath = fullfile(defaultFolder,'!Summary');
mkdir(savePath);

figurePath = fullfile(savePath,'Summary Figures');
mkdir(figurePath);

saveFigures(figurePath,true,0);