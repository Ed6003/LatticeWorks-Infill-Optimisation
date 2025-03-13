% Clear environment and set default settings
clear; close all; clc; 

% Define plot settings
setPlotSettings();

% Initialise figures and parameters
count = 1;
nCols = 3; nRows = 2; linspace_steps = 64;
figures = createFigures(10);  % Create 10 figures
original_points = linspace(1.3, -1.3, linspace_steps);
selected_idx = round(linspace(1, linspace_steps, nCols * nRows));

% Process selected density values
for density in original_points[selected_idx]
    resultStruct = loadSimulationResults(density);

    % Generate figures for isosurface and displacement
    plotIsoSurface(figures(1), resultStruct);
    plotOrthographicProjection(figures(2), resultStruct);
    plotIsoSurfaceGeneration(figures(3), resultStruct);
    plotGeogramRemeshed(figures(4), resultStruct);
    plotNodeSets(figures(5), resultStruct);
    plotBoundaryConditions(figures(6), resultStruct);
    plotTrackedDisplacement(figures(7), resultStruct.xz, resultStruct.Uz, 'Z');
    plotTrackedDisplacement(figures(8), resultStruct.xy, resultStruct.Uy, 'Y');
    plotPoissonsRatio(figures(9), resultStruct);
    plotYoungsModulus(figures(10), resultStruct);

    count++;
end

% Save summary figures
saveSummaryFigures('D:\TechnicalReport\LatticeProperties_LinearFinal\!Summary');

% Function Definitions
function setPlotSettings()
    global figureTitles; 
    figureTitles = true;
    % Define colour map and styles
end

function figures = createFigures(num)
    for i = 1:num
        figures(i) = cFigure();
    end
end

function resultStruct = loadSimulationResults(density)
    folderPath = 'D:\TechnicalReport\LatticeProperties_LinearFinal';
    savePath = fullfile(folderPath, sprintf('%.5g_Lattice_Density', density));
    matlabPath = fullfile(savePath, 'simulation_results.mat');
    load(matlabPath);
end

function plotIsoSurface(fig, resultStruct)
    figure(fig);
    subplot(nRows, nCols, count);
    title(sprintf('Isovalue: %.3g', resultStruct.level_set));
    hp1 = gpatch(resultStruct.f, resultStruct.v, resultStruct.c, 'none', 1);
    hp1.FaceColor = 'flat';
    colormap(gca, gjet(6));
    axisGeom(gca, fontSize); camlight headlight;
end

function plotOrthographicProjection(fig, resultStruct)
    figure(fig);
    subplot(nRows, nCols, count);
    title(sprintf('Isovalue: %.3g', resultStruct.level_set));
    hp2 = gpatch(resultStruct.f, resultStruct.v, 'none', 1);
    hp2.FaceColor = [0 0.4470 0.7410];
    axis equal; xlim([min(resultStruct.v(:,1)), max(resultStruct.v(:,1))]);
    ylim([min(resultStruct.v(:,2)), max(resultStruct.v(:,2))]);
end

function plotTrackedDisplacement(fig, xData, yData, direction)
    figure(fig);
    subplot(nRows, nCols, count);
    title(sprintf('Infill: %.3g%%', resultStruct.infill_percentage));
    scatter(xData, yData, 'x', 'MarkerEdgeColor', direction == 'Z' ? blue : cyan);
    xlabel('Global X Coordinate (mm)');
    ylabel(sprintf('%s Displacement (mm)', direction));
    grid on;
end

function plotPoissonsRatio(fig, resultStruct)
    figure(fig);
    subplot(nRows, nCols, count);
    title(sprintf('Infill: %.3g%%', resultStruct.infill_percentage));
    plot(resultStruct.U_x, resultStruct.poisson.poisson_xz_mean, '-', 'Color', blue);
    plot(resultStruct.U_x, resultStruct.poisson.poisson_xy_mean, '-', 'Color', cyan);
    xlabel('Applied X Displacement (mm)');
    ylabel('Poisson Ratio');
end

function plotYoungsModulus(fig, resultStruct)
    figure(fig);
    subplot(nRows, nCols, count);
    title(sprintf('Infill: %.3g%%', resultStruct.infill_percentage));
    plot(abs(resultStruct.strain.strain_y_mean), resultStruct.stress, '-', 'Color', blue);
    plot(abs(resultStruct.strain.strain_z_mean), resultStruct.stress, '-', 'Color', cyan);
    xlabel('Strain');
    ylabel('Stress (MPa)');
end

function saveSummaryFigures(folderPath)
    mkdir(folderPath);
    saveFigures(folderPath, true, 0);
end