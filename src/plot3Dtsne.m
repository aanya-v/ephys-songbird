function plot3Dtsne(tsneResult,daysLabels,sylsLabels)
%UNTITLED2 Summary of this function goes here

% Assign unique colors based on the first label set
    uniqueLabels = unique(daysLabels);
    numUniqueLabels = numel(uniqueLabels);
    colors = lines(numUniqueLabels); % You can use any colormap you prefer
    
    % Create a color map based on labels
    labelColors = zeros(length(daysLabels), 3);
    for i = 1:numUniqueLabels
        labelColors(daysLabels == uniqueLabels(i), :) = repmat(colors(i, :), [sum(daysLabels == uniqueLabels(i)), 1]);
    end
    
% % Step 5b: Optional - assign unique colors based on second label set
    uniqueLabels2 = unique(sylsLabels);
    numUniqueLabels2 = numel(uniqueLabels2);
    colors2 = lines(numUniqueLabels2); % You can use any colormap you prefer

    % Create a color map based on labels
    labelColors2 = zeros(length(sylsLabels), 3);
    for i = 1:numUniqueLabels2
        labelColors2(sylsLabels == uniqueLabels2(i), :) = repmat(colors2(i, :), [sum(sylsLabels == uniqueLabels2(i)), 1]);
    end

%% Create 3D scatter plot
figure;
hold on;
    for i = 1:numUniqueLabels
        idx = daysLabels == uniqueLabels(i);
        scatter3(tsneResult(idx, 1), tsneResult(idx, 2), tsneResult(idx, 3), 3, labelColors(idx, :), 'filled');
        
        % Label each cluster with the associated character
        centroid = mean(tsneResult(idx, :), 1);
    end

for i = 1:numUniqueLabels2
    idx = sylsLabels == uniqueLabels2(i);
    % scatter3(tsneResult(idx, 1), tsneResult(idx, 2), tsneResult(idx, 3), 3, labelColors2(idx, :), 'filled');
    
    % Label each cluster with the associated character
    centroid2 = mean(tsneResult(idx, :), 1);
    text(centroid2(1), centroid2(2), centroid2(3), uniqueLabels2(i), 'FontSize', 30, 'FontWeight', 'bold');
end
hold off;

% title(['br177yw112 - syllables: ' syls ' by day']);
% subtitle('settings: exact, euclidean')
xlabel('t-SNE1');
ylabel('t-SNE2');
zlabel('t-SNE3');

% Create string array for the legend
dayStrings = arrayfun(@(x) ['D' num2str(x)], uniqueLabels, 'UniformOutput', false);
legend(dayStrings, "Location","southeast")
view(3)
set(gcf,'Visible','On')

end