function exportNamedFigures(numbers, figureNames, projectDirectory)
%EXPORTNAMEDFIGURES  Export named MATLAB figures to plots/ or SGN_OUTPUT_DIR.
%
%   Writes each requested figure number to its configured PNG file name
%   under the project plots folder, or under SGN_OUTPUT_DIR when set.
%
%   Author: Pasquale Marzaioli

outputDirectory = getenv('SGN_OUTPUT_DIR');
if isempty(outputDirectory)
    outputDirectory = fullfile(projectDirectory, 'plots');
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
for number = numbers
    % Fail early if a new figure has no meaningful external file name.
    assert(number <= numel(figureNames) && ~isempty(figureNames{number}), ...
        'Missing export file name for figure %d.', number);
    drawnow;
    exportgraphics(figure(number), fullfile(outputDirectory, ...
        figureNames{number}), 'Resolution', 160);
end
end
