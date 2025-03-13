% Initialise environment
clear all variables; 
close all figures; 
clc; 

% Settings
updateSim = true;              % Flag to recompute Abaqus simulation
overwrite = false;             % Flag to overwrite existing simulated files
coarseMesh = false;            % Flag for coarse mesh testing
figureTitles = true;           % Flag to show figure titles
loadSimDataFromMat = false;    % Flag to load data from .mat file if sim exists
doSaveFigures = true;          % Flag to save figures

% Change directory to the simulation path
cd('C:\Users\rusco\OneDrive - University of Warwick\Admin\Archives\Documents\GitHub\LatticeWorks');

% Start CPU logging
system('logman start CPU_Log');

% Define plot settings
cMap = generateColourMap();     % Custom function to generate colour map
fontSize = 15;
markerSize = 20;

% Control parameters
pointSpacing = 0.15;            % Point spacing for mesh generation

% Geometry generation loop
for i = linspace(1.3,-1.3,64) % Loop over isosurface levels

    nFigures = 0;                % Initialize figure count
    levelset = i;               % Isosurface level

    % Geometry dimensions
    x_length = 30; 
    y_length = 10; 
    z_length = 10; 

    % Input structure for mesh generation
    inputStruct = initialiseInputStruct(x_length, y_length, z_length);

    % Prepare save path for results
    savePath = generateSavePath(i);

    % Check if simulation is already completed
    if simulationCompleted(savePath) && ~overwrite
        fprintf('[%.5g] Lattice density already completed, skipping ...\n', i);
        continue;  % Skip to next iteration
    end

    % Ensure output directory exists
    createDirectory(savePath);

    % Prepare Abaqus input file names
    abaqusInpFileName = createAbaqusInputFileName(savePath);
    abaqusDATFileName = createAbaqusDatFileName(savePath);

    % Material properties
    resultStruct.material = setMaterialProperties();

    % Generate geometry
    [S,X,Y,Z] = gradTPMS(inputStruct);
    [F,V] = isosurface(X,Y,Z,S,levelset);
    [fc,vc] = isocaps(X,Y,Z,S,levelset, 'above');
    [f,v,c] = FV_arrange(F,V,fc,vc);  % Join and clean geometry

    % Visualise surface
    visualiseSurface(f,v,c,figureTitles,fontSize);

    % Remesh geometry
    optionStruct.pointSpacing = pointSpacing;
    [F,V] = ggremesh(f,v,optionStruct);
    resultStruct.F = F; 
    resultStruct.V = V;

    % Tetrahedral meshing
    performTetGenMeshing(inputStruct, F, V, coarseMesh, optionStruct);

    % Visualise mesh
    meshView(meshOutput);

    % Node labels and boundary condition selection
    C_vertex = labelBoundaryNodes(V, x_length, y_length, z_length);
    resultStruct.C_vertex = C_vertex;

    % Visualise boundary conditions
    visualiseBoundaryConditions(Fb, V, C_vertex, markerSize, figureTitles, fontSize);

    % Save all figures
    saveFigures(savePath, doSaveFigures, nFigures);

    % Prepare Abaqus input structure
    abaqus_spec = prepareAbaqusInputStruct(inputStruct, V, E_youngs, v_poisson, bcEncastreList, bcLoadList);

    % Write Abaqus input file
    if updateSim
        abaqusStruct2inp(abaqus_spec, abaqusInpFileName);
        runAbaqusJob(abaqusInpFileNamePart, savePath);
    end

    % Import and visualise Abaqus results
    abaqusData = importAbaqusResults(loadSimDataFromMat, saveFile, abaqusDATFileName);

    % Fetch element data from results
    [E_effectiveStress, E_effectiveStrain] = fetchElementData(abaqusData, E);

    % Animate deformations
    animateDeformations(abaqusData, V, fontSize);

    % Calculate and visualise Poisson's ratio
    calculatePoissonsRatio(abaqusData, resultStruct);

    % Finalize results and save
    resultStruct = finalizeResults(resultStruct, abaqusData, E_effectiveStress, E_effectiveStrain);
    saveSimulationResults(savePath, resultStruct);

    % Stop CPU logging
    system('logman stop CPU_Log');
end