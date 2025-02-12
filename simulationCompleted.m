function isCompleted = simulationCompleted(parentDir)
% simulationCompleted  Check if a given simulation directory has been completed.
%   isCompleted = simulationCompleted(parentDir) returns true if the base
%   directory (parentDir) contains all the required simulation files and
%   if the Figures subfolder contains all the required figure files.

    % Define the required simulation files (expected to be in parentDir)
    requiredSimulationFiles = {
        'Lattice_FEA.com', 'Lattice_FEA.dat', 'Lattice_FEA.inp', 'Lattice_FEA.msg', ...
        'Lattice_FEA.odb', 'Lattice_FEA.prt', 'Lattice_FEA.sta', 'simulation_results.mat'
    };

    % Define the required figure files (expected to be in parentDir\Figures)
    requiredFigureFiles = {
        'Figure_1.fig', 'Figure_1.png', 'Figure_2.fig', 'Figure_2.png', ...
        'Figure_3.fig', 'Figure_3.png', 'Figure_4.fig', 'Figure_4.png', ...
        'Figure_5.fig', 'Figure_5.png', 'Figure_6.fig', 'Figure_6.png', ...
        'Figure_7.fig', 'Figure_7.png', 'Figure_8.fig', 'Figure_8.png', ...
        'Figure_9.fig', 'Figure_9.png', 'Figure_10.fig', 'Figure_10.png', ...
        'Figure_11.fig', 'Figure_11.png', 'Figure_12.fig', 'Figure_12.png', ...
        'Figure_13.fig', 'Figure_13.png', 'Figure_14.fig', 'Figure_14.png', ...
        'Figure_15.fig', 'Figure_15.png', ...
    };

    % Check if the base directory exists
    if ~isfolder(parentDir)
        isCompleted = 0;
        return;
    end

    % Check if all required simulation files exist in the base directory.
    simulationFilesExist = all(cellfun(@(f) exist(fullfile(parentDir, f), 'file') == 2, requiredSimulationFiles));

    % Define the Figures subfolder path and check if it exists.
    figuresDir = fullfile(parentDir, 'Figures');
    if ~isfolder(figuresDir)
        figureFilesExist = false;
    else
        % Check if all required figure files exist in the Figures subfolder.
        figureFilesExist = all(cellfun(@(f) exist(fullfile(figuresDir, f), 'file') == 2, requiredFigureFiles));
    end

    isCompleted = simulationFilesExist && figureFilesExist;
end