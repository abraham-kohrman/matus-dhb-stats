function [fake1, fake2] = randomizePoints(real1, real2)
    %get sizes    
    len1 = length(real1);
    len2 = length(real2);
    fullen = len1 + len2;
    
    %concatenate
    stacked = horzcat(real1, real2);
    
    %randomize 
    x = randperm(fullen);
    randidx1 = x(1:len1);
    randidx2 = x(len1+1:end);
    
    %return
    fake1 = stacked(randidx1);
    fake2 = stacked(randidx2);
   
end