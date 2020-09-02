function areaOut = deltaFunction(curve1, curve2)
%a very basic delta function which simply subtracts the sum of each curve. 
    areaOut = nansum(curve1(:)) - nansum(curve2(:)); %this is a lame way of implementing this--but it's a start!
    
end
