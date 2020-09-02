function [randomized1, randomized2] = randomizeTraces_new(full_set,len1,replace)
    %
    ncol = size(full_set,2); 
    if replace
        x= randi(ncol,1,ncol);
    else
        x = randperm(ncol);
    end
    randidx1 = x(1:len1);
    randidx2 = x(len1+1:end);
    
    randomized1 = full_set(:,randidx1);
    randomized2 = full_set(:,randidx2);
end