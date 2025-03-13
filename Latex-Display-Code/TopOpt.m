% Set working directory and save path
cd("C:\Users\rusco\OneDrive - University of Warwick\Admin\Archives\Documents\GitHub\LatticeWorks\TopOpt");
savePath = 'D:\TechnicalReport\Variable-TPMS-Figures\TopOpt';

% Load model data
load('sample.mat');
model = createpde("structural", "static-solid");
importGeometry(model, 'Planet-18-Teeth-1.8mm-Mod.step');
scale(model.Geometry, 1e3);  % Scale to mm

% Set material properties
setMaterialProperties(model, 50e6, 1666e6, 0.38);  % (Yield strength, E, v)

% Visualise geometry
visualiseGeometry(model);

% Define parameters for teeth and load
teeth_faces = [14,24,31,38,48,55,62,72,79,86,96,103,110,120,127,134,144,151];
pressure = calculatePressure(10, 15e-3, teeth_faces);

% Apply loads and constraints
applyBoundaryConditions(model, teeth_faces, pressure);

% Generate mesh and solve
generateMesh(model, 'GeometricOrder', 'quadratic', 'Hmax', 0.0005);
result = solve(model);

% Display results
displayResults(result);

% Interpolate stress data
[xq, yq, zq, stressGrid] = interpolateStressData(result);

% Mask and visualise stress distribution
maskAndVisualiseStress(xq, yq, zq, stressGrid);

% Create and visualise TPMS
createAndVisualiseTPMS(xq, yq, zq, stressGrid, savePath);

% Save figures
saveFigures(savePath, true, 0);

% Function Definitions
function setMaterialProperties(model, yield_strength, E, v)
    % Sets the material properties for the model
    structuralProperties(model, "YoungsModulus", E, "PoissonsRatio", v);
end

function visualiseGeometry(model)
    cFigure();
    pdegplot(model, 'FaceLabels', 'on', 'FaceAlpha', 0.5);
    title('Imported Gear Geometry');
end

function pressure = calculatePressure(torque, teeth_distance, teeth_faces)
    teeth_force = (torque / teeth_distance) / numel(teeth_faces);  % N
    area_mm2 = 66.441 * 1e-6;  % mm^2
    pressure = teeth_force / area_mm2;  % Pa
end

function applyBoundaryConditions(model, teeth_faces, pressure)
    structuralBoundaryLoad(model, "Face", teeth_faces, "Pressure", pressure);
    structuralBC(model, "Face", 8, "Constraint", "fixed");
end

function displayResults(result)
    fprintf('Max Von Mises: %g MPa.\n', max(result.VonMisesStress) / 1e6);
end

function [xq, yq, zq, stressGrid] = interpolateStressData(result)
    % Interpolates Von Mises stress over a grid
    x = result.Mesh.Nodes(1, :);
    y = result.Mesh.Nodes(2, :);
    z = result.Mesh.Nodes(3, :);
    vonMisesStress = result.VonMisesStress;

    n = 200;
    [Xq, Yq, Zq] = meshgrid(linspace(min(x), max(x), n), linspace(min(y), max(y), n), linspace(min(z), max(z), round((max(z) - min(z))/(max(x)-min(x))/(n-1)) + 1));
    F = scatteredInterpolant(x', y', z', vonMisesStress, 'linear', 'none');
    stressGrid = F(Xq, Yq, Zq);
    stressGrid(isnan(stressGrid)) = 0;
end

function maskAndVisualiseStress(xq, yq, zq, stressGrid)
    % Apply masking and visualisation of interpolated stress
    alphaVal = 1;  % Threshold
    shp = alphaShape(xq', yq', z', alphaVal);
    insideMask = inShape(shp, xq, yq, zq);
    stressGrid(~insideMask) = NaN;

    % Create surface plot
    cFigure();
    surf(xq(:,:,1), yq(:,:,1), stressGrid(:,:,1), 'EdgeColor', 'none');
    colorbar();
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    zlabel('Von Mises Stress');
end

function createAndVisualiseTPMS(xq, yq, zq, stressGrid, savePath)
    % Generate and visualise TPMS structure
    % Logic for TPMS generation and visualisation goes here

    % Save generated mesh
    TR = triangulation(Fsn, Vsn);
    stlwrite(TR, fullfile(savePath, 'OutputGear.stl'));
end