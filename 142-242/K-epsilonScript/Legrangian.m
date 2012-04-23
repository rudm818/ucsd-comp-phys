function [L] = Legrangian(idxI,idxJ)
%LEGRANGIAN Summary of this function goes here
%   Generates Legrangian elements of K-epsilon matrix
    M = 1;
    Domain = 10;
    Range  = 10;
    lengthUnit = 1; %scaling
    
    
    % map index space to x space
    slope = (2*Range)/(Domain-1);
    intercept = -(Domain+1);
    XposI = lengthUnit*(slope*idxI+intercept);
    XposJ = lengthUnit*(slope*idxJ+intercept);
    
    L = M/2*(XposI-XposJ)^2 - 1/2*((XposI+XposJ)/2)^2; % T(x)-V(x)

end

