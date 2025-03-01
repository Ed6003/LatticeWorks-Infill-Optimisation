%(*@\draftcomment{the prefix for every sub-heading in this script is "savefigures"}@*)
function nFigures = saveFigures(savePath,doSaveFigures,nFigures)
    figsavePath = fullfile(savePath, 'Figures');
    figHandles = findobj('Type', 'figure');

    if ~isfolder(figsavePath)
        mkdir(figsavePath);
    end

    %(*@\codesubsection{Iteratively Save Figures}{savefigures-iteratively-save-figures}@*)
    if doSaveFigures
        for j = 1:length(figHandles)

            figHandle = figHandles(j);
            
            % generate unique filename
            figName = sprintf('Figure_%d', figHandle.Number + nFigures);
            
            savefig(figHandle, fullfile(figsavePath, [figName, '.fig'])); % save as .fig  
            saveas(figHandle, fullfile(figsavePath, [figName, '.png'])); % save as .png
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