function [num_time_bins] = calcNumBins(spectrograms,syl_duration)

if syl_duration =="longest"
    maxColumns = 0;
    for i = 1:length(spectrograms)
        for j = 1:length(spectrograms{i})
            currentColumns = size(spectrograms{i}, 2);
            if currentColumns > maxColumns
                maxColumns = currentColumns;
            end
        end
    end
    
    num_time_bins = maxColumns;

elseif syl_duration == "shortest"
    minColumns = 1000000;
    for i = 1:length(spectrograms)
        for j = 1:length(spectrograms{i})
            currentColumns = size(spectrograms{i}, 2);
            if currentColumns < minColumns
                minColumns = currentColumns;
            end
        end
    end
    
    num_time_bins = minColumns;
end
end