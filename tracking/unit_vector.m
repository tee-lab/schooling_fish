function [unitvectrs] = unit_vector(pvel)
univec = zeros(size(pvel,1),2);
    
        for j=1:size(pvel,1)
            if norm(pvel(j,:)) ~= 0
                univec(j,:) = pvel(j,:)./norm(pvel(j,:)) ; 
            else
                univec(j,:) = [0 0];
            end
        
        end
        unitvectrs = univec ; 
end