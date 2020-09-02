function [mean1,mean2] = fakeMeans(key1, key2)
    %pull the relevant datasets
    
    temp_mean1 = key1;
    temp_mean2 = key2;
    temp_mean1(temp_mean1 <= 0 ) = NaN;
    temp_mean2(temp_mean2 <= 0 ) = NaN; 
    
    %get their lengths
    [len1,~] = size(temp_mean1);
    [len2,~] = size(temp_mean2);
    
    %figure out which is bigger
    if len1>len2
        len = len2;
    else
        len = len1;
    end
    
    %clip to shorer length--> replace with interpolation!
    temp_mean1 = temp_mean1(1:len,2:end);
    temp_mean2 = temp_mean2(1:len,2:end);
    
    %get the mean in the second dimension
    mean1 = nanmean(temp_mean1,2);
    mean2 = nanmean(temp_mean2,2);
    
end
