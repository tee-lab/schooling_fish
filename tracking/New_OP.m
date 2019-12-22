function [OP1,OP2,velx,vely] = New_OP(data_sheet,ini,maxframe,count,p)
%%
ordr_prmtr1= zeros(maxframe,1);
velx= zeros(maxframe,1);
vely= zeros(maxframe,1);
% c=0;
parfor i=ini:maxframe
%     c=c+1;
    ptest= data_sheet(data_sheet(:,1) == i, :); %get positions of ith frame
    
    %% Filter lesser nbr of detected indi. 
    
%     [nz,~] =find(~isnan(ptest(:,7)) & ~isnan(ptest(:,8)));
%     if nbr_of_indi == 15d
%     percent = length(nz)/nbr_of_indi;
%     else
%         percent = 100;
%     end
    %%
    if size(ptest,1) > 0 %&& count(i,3) > p
    
        pos = ptest(~isnan(ptest(:,7)),3:4); % get only positions
        pos1 = pos(:,1);
        pos2 = pos(:,2);
    % For Center of mass 
        if size(pos,1) > 1
            centr_mas = mean(pos);
        else
            centr_mas = pos ;
        end
    
    
    pos_vec1 = pos1 - centr_mas(:,1);
    pos_vec2 = pos2 - centr_mas(:,2);
    pos_vec = cat(2,pos_vec1,pos_vec2);
    %%
    Pvel = ptest(:,7:8);
    pvel=Pvel(~isnan(Pvel(:,1)),:);
    univec = unit_vector(pvel);
    pos_univec = unit_vector(pos_vec);

%         ordr_prmtr(i,2) = sqrt( (sum(univec(:,1))).^2 + (sum(univec(:,2))).^2 )  / length(univec);
        velx(i) = sum(univec(:,1)) / length(univec);
        vely(i) = sum(univec(:,2)) / length(univec);
        ordr_prmtr1(i,1) =  norm(sum(univec))/ length(univec);
        ordr_prmtr2(i,1) =  norm(sum(univec.*pos_univec))/ length(univec);
%         ordr_prmtr2(i,1) = norm(sum(pvel)/length(pvel));
%         ordr_prmtr3(c,1:2) = centr_mas ;
        i
    end
end
OP1 = ordr_prmtr1;
OP2 = ordr_prmtr2;
end
