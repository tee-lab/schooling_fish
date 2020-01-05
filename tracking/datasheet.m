%for saving in order
function [data_sheet] = datasheet(Q_loc_estimateX,Q_loc_estimateY,steps,frame_rate)
clear data_sheet;
% nbr_of_animals=100;
% data_sheet=nan(nbr_of_animals*length(Q_loc_estimateX(:,1)),8);
% To make a vector components from the tracked data:
clear Vx Vy;

nbrof_frames=length(Q_loc_estimateX(:,1));
% parfor i=1:floor(nbrof_frames-1)
% Vx(i,:)=(Q_loc_estimateX(i+1,:)-Q_loc_estimateX(i,:))/(steps/30);
% Vy(i,:)=(Q_loc_estimateY(i+1,:)-Q_loc_estimateY(i,:))/(steps/30);
% end
%------------------------------------------
Q_loc_estimateX1=Q_loc_estimateX(2:end,:);
Q_loc_estimateY1=Q_loc_estimateY(2:end,:);

Vx=(Q_loc_estimateX1-Q_loc_estimateX(1:length(Q_loc_estimateX1(:,1)),:))/(steps/frame_rate);
Vy=(Q_loc_estimateY1-Q_loc_estimateY(1:length(Q_loc_estimateX1(:,1)),:))/(steps/frame_rate);

clear Q_loc_estimateX1 Q_loc_estimateY1 ;
%-------------------------------------------

clear vx vy;
vx(2:length(Vx(:,1))+1,1:length(Vx(1,:)))=Vx;
vy(2:length(Vx(:,1))+1,1:length(Vx(1,:)))=Vy;
clear Vx Vy;
Vx=vx; Vy=vy;

d=0;

%% sort wrt frame nbr

for i=1:length(Vx(:,1))
    c=0;
    for j=1:length(Vx(1,:))
            if (Q_loc_estimateX(i,j)) > 0
            d=d+1;
            data_sheet(d,1)=i; %frame nbr index, you need to multiply by steps to get actual frame number. 
                if (Q_loc_estimateX(i,j)) > 0
                c=c+1;
                data_sheet(d,2)=c; %identity
                data_sheet(d,3)=Q_loc_estimateX(i,j); % x position
                data_sheet(d,4)=Q_loc_estimateY(i,j); % y position
                    
                    if i==1
                        data_sheet(d,5)=0;
                        data_sheet(d,6)=0;
                    else
                        data_sheet(d,5)=Q_loc_estimateX(i-1,j); % x position previous
                        data_sheet(d,6)=Q_loc_estimateY(i-1,j); % y position previous
                    end
                data_sheet(d,7)=Vx(i,j); %velocity x-component
                data_sheet(d,8)=Vy(i,j); %velocity y-component
                else
                c=c;
                end
            end
    end
     
end

end
% avg_velocity;
% NewOP;
% NewOP2;