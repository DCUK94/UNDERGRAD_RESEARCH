function [ out ] = checkForPixelEdges( nhood )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 if(nhood(2,2) == 1)
    sum_4_nhood = nhood(1,2) + nhood(2,1) + nhood(3,2) + nhood(2,3);
  
    if (sum_4_nhood == 3)
        out = 1; 
    elseif (sum_4_nhood == 2)
        out = 2;
    elseif (sum_4_nhood == 1)
        out = 3;
    elseif((sum(nhood(:))-nhood(2,2)) == 0)
        out = 4;
    else
         out = 0;
    end
else    
    out = 0;
end
end

