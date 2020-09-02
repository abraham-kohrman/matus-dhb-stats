function areaOut = deltaFunctionLists(curve1, curve2)
%this is a new deta function which does pairwise subtraction
%     temp = curve1 - curve2;
%     temp = abs(temp);
%     areaOut = nansum(temp);
    areaOut = nansum(abs(curve1 - curve2));
end