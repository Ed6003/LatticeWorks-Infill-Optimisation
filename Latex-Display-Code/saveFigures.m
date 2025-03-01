function nFigures = saveFigures(savePath,doSaveFigures,nFigures)
    figsavePath = fullfile(savePath, 'Figures');
    figHandles = findobj('Type', 'figure');

    if ~isfolder(figsavePath)
        mkdir(figsavePath);
    end

    
    if doSaveFigures
        % Iterate over each figure
        for j = 1:length(figHandles)
            % Get the current figure handle
            figHandle = figHandles(j);
            
            % Generate a unique filename based on the figure number
            figName = sprintf('Figure_%d', figHandle.Number + nFigures);
            
            % Save the figure as a .fig file
            savefig(figHandle, fullfile(figsavePath, [figName, '.fig']));
            
            % Save the figure as a .png file
            saveas(figHandle, fullfile(figsavePath, [figName, '.png']));
            
            % Close the figure
            close(figHandle);
        end
    else
        for j = 1:length(figHandles)
            figHandle = figHandles(j);
            close(figHandle);
        end
    end

    nFigures = length(figHandles);
end