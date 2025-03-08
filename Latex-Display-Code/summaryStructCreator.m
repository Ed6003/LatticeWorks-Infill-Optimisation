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

%(*@
for i = 1:numel(summaryStruct)
    summaryStruct(i).youngs_modulus.youngs_modulus_xz_mean = summaryStruct(i).youngs_modulus.youngs_modulus_xz_mean / 10;
    summaryStruct(i).youngs_modulus.youngs_modulus_xz_median = summaryStruct(i).youngs_modulus.youngs_modulus_xz_median / 10;
    summaryStruct(i).youngs_modulus.youngs_modulus_xy_mean = summaryStruct(i).youngs_modulus.youngs_modulus_xy_mean / 10;
    summaryStruct(i).youngs_modulus.youngs_modulus_xy_median = summaryStruct(i).youngs_modulus.youngs_modulus_xy_median / 10;
end
%@*)

%(*@\codesubsection{Figure 1}{summarystruct-figure-1}@*)
cFigure;
hold on; grid on;
plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xy_mean),'Color',blue, 'LineWidth', 1.5)
plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xy_median), '--','Color',blue, 'LineWidth', 1.5)
plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xz_mean),'Color',cyan, 'LineWidth', 1.5)
plot(vertcat(summaryStruct.infill_percentage), vertcat(summaryStruct.poisson_xz_median), '--','Color',cyan, 'LineWidth', 1.5)


title("Calculated Poisson's Ratio",'FontSize',14)
legend('XY Mean', 'XY Median', 'XZ Mean', 'XZ Median', 'Location','best');
xlabel("Infill Percentage (%)")
ylabel("Poisson's Ratio")

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
% preallocate
n = numel(summaryStruct);
infill = zeros(n,1);
poisson_xy_median_val = zeros(n,1);
poisson_xy_mad = zeros(n,1);
poisson_xz_median_val = zeros(n,1);
poisson_xz_mad = zeros(n,1);

for i = 1:n
    infill(i) = summaryStruct(i).infill_percentage;
    dataXy = [summaryStruct(i).poisson_xy_mean, summaryStruct(i).poisson_xy_median];
    poisson_xy_median_val(i) = median(dataXy);
    poisson_xy_mad(i) = mad(dataXy, 1);
    dataXz = [summaryStruct(i).poisson_xz_mean, summaryStruct(i).poisson_xz_median];
    poisson_xz_median_val(i) = median(dataXz);
    poisson_xz_mad(i) = mad(dataXz, 1);
end

figure;
hold on; grid on;

hXy = errorbar(infill, poisson_xy_median_val, poisson_xy_mad, 'o-', ...
    'LineWidth', 1.5, 'CapSize', 5, 'Color', blue);

% create 95% heuristic region for xy
xFill = [infill; flipud(infill)];
ciXyUpper = poisson_xy_median_val + 1.96 * poisson_xy_mad; % equal to 95% confidence
ciXyLower = poisson_xy_median_val - 1.96 * poisson_xy_mad;
yFillXy = [ciXyLower; flipud(ciXyUpper)];
fill(xFill, yFillXy, blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot XZ with error bars using MAD (cyan)
hXz = errorbar(infill, poisson_xz_median_val, poisson_xz_mad, 'o-', ...
    'LineWidth', 1.5, 'CapSize', 5, 'Color', cyan);

% create xz heuristic
ciXzUpper = poisson_xz_median_val + 1.96 * poisson_xz_mad;
ciXzLower = poisson_xz_median_val - 1.96 * poisson_xz_mad;
yFillXz = [ciXzLower; flipud(ciXzUpper)];
fill(xFill, yFillXz, cyan, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

title("Calculated Poisson's Ratio with Robust Uncertainty", 'FontSize', 14);
xlabel("Infill Percentage (%)");
ylabel("Poisson's Ratio");
legend([hXy, hXz], 'XY Median ± MAD', 'XZ Median ± MAD', 'Location', 'best');
hold off;
%%
%(*@\codesubsection{Figure 11}{summarystruct-figure-11}@*)
youngsModulus = vertcat(summaryStruct.youngs_modulus);
infillPercentage = vertcat(summaryStruct.infill_percentage);

n = numel(summaryStruct);
youngs_modulus_xy_median = zeros(n, 1);
youngs_modulus_xy_mad = zeros(n, 1);
youngs_modulus_xz_median = zeros(n, 1);
youngs_modulus_xz_mad = zeros(n, 1);

for i = 1:n
    dataXy = [youngsModulus(i).youngs_modulus_xy_mean, youngsModulus(i).youngs_modulus_xy_median];
    youngs_modulus_xy_median(i) = median(dataXy);
    youngs_modulus_xy_mad(i) = mad(dataXy, 1);
    
    dataXz = [youngsModulus(i).youngs_modulus_xz_mean, youngsModulus(i).youngs_modulus_xz_median];
    youngs_modulus_xz_median(i) = median(dataXz);
    youngs_modulus_xz_mad(i) = mad(dataXz, 1);
end

figure;
hold on; grid on;

hXy = errorbar(infillPercentage, youngs_modulus_xy_median, youngs_modulus_xy_mad, 'o-', ...
    'LineWidth', 1.5, 'CapSize', 5, 'Color', blue);
hXz = errorbar(infillPercentage, youngs_modulus_xz_median, youngs_modulus_xz_mad, 'o-', ...
    'LineWidth', 1.5, 'CapSize', 5, 'Color', cyan);

xFill = [infillPercentage; flipud(infillPercentage)];

ciXyLower = youngs_modulus_xy_median - 1.96 * youngs_modulus_xy_mad;
ciXyUpper = youngs_modulus_xy_median + 1.96 * youngs_modulus_xy_mad;
yFillXy = [ciXyLower; flipud(ciXyUpper)];
fill(xFill, yFillXy, blue, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

ciXzLower = youngs_modulus_xz_median - 1.96 * youngs_modulus_xz_mad;
ciXzUpper = youngs_modulus_xz_median + 1.96 * youngs_modulus_xz_mad;
yFillXz = [ciXzLower; flipud(ciXzUpper)];
fill(xFill, yFillXz, cyan, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

title("Calculated Young's Moduli with Uncertainty", 'FontSize', 14);
legend([hXy, hXz], 'XY Median ± MAD', 'XZ Median ± MAD', 'Location', 'best');
xlabel("Infill Percentage (%)");
ylabel("Young's Modulus (MPa)");
hold off;

%%
savePath = fullfile(defaultFolder,'!Summary');
mkdir(savePath);

figurePath = fullfile(savePath,'Summary Figures');
mkdir(figurePath);

saveFigures(figurePath,true,0);