function paramOut = deltaFunction_points(curve1, curve2,type)
    switch type
        case 'Mean'
            paramOut = nanmean(curve1) - nanmean(curve2);
        case 'Median'
            paramOut = nanmedian(curve1) - nanmedian(curve2);
        case 'First'
            paramOut = curve1(1) - curve2(1);
        otherwise
            paramOut = -1;
    end
end