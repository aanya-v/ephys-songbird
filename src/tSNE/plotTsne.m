function plotTsne(tsneResult, dimension, daysLabels, sylsLabels, colorBy)
    % Check if the user has chosen to color by day or syllable labels
    if nargin < 5
        colorBy = 'day'; % Default color by day
    end

    % Assign unique colors based on the chosen label set
    if strcmp(colorBy, 'day')
        labelsToColor = daysLabels;
        % uniqueLabels = unique(daysLabels);
        % labelColors = assignColors(daysLabels, uniqueLabels);
    else
        labelsToColor = sylsLabels;
        % uniqueLabels = unique(sylsLabels);
        % labelColors = assignColors(sylsLabels, uniqueLabels);
    end

    uniqueLabels = unique(labelsToColor);
    labelColors = assignColors(labelsToColor, uniqueLabels);

    % Create scatter plot based on dimension
    figure;
    hold on;

    if dimension == 3
        % 3D scatter plot
        scatterPlot(tsneResult, labelColors, dimension, uniqueLabels, labelsToColor,sylsLabels, colorBy);
        zlabel('t-SNE3');
        view(3);
    else
        % 2D scatter plot
        scatterPlot(tsneResult, labelColors, dimension, uniqueLabels, labelsToColor, sylsLabels, colorBy);
    end

    % Common plot settings
    xlabel('t-SNE1');
    ylabel('t-SNE2');
    legend(createLegendStrings(uniqueLabels), "Location", "southeast");
    hold off;
    set(gcf, 'Visible', 'On');
end

function labelColors = assignColors(labels, uniqueLabels)
    numUniqueLabels = numel(uniqueLabels);
    colors = turbo(numUniqueLabels);
    labelColors = zeros(length(labels), 3);
    for i = 1:numUniqueLabels
        labelColors(labels == uniqueLabels(i), :) = repmat(colors(i, :), sum(labels == uniqueLabels(i)), 1);
    end
end

function scatterPlot(tsneResult, labelColors, dimension, uniqueLabels, labelsToColor, sylsLabels, colorBy)
    for i = 1:numel(uniqueLabels)
        
        idx = (labelsToColor == uniqueLabels(i));
        if dimension == 3
            scatter3(tsneResult(idx, 1), tsneResult(idx, 2), tsneResult(idx, 3), 12, labelColors(idx, :), 'filled');
        else
            scatter(tsneResult(idx, 1), tsneResult(idx, 2), 12, labelColors(idx, :), 'filled');
        end
        % Add centroids for the other label type
        if strcmp(colorBy, 'day') && dimension == 3
            centroid = mean(tsneResult(sylsLabels == uniqueLabels(i), :), 1);
            text(centroid(1), centroid(2), centroid(3), uniqueLabels(i), 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
end

function legendStrings = createLegendStrings(uniqueLabels)
    legendStrings = arrayfun(@(x) ['Label ' num2str(x)], uniqueLabels, 'UniformOutput', false);
end