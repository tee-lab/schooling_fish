%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jitesh Jhawar                                               %
% Description:                                                        %
%                                                                     %
% This code reads through meta_data and generates the order parameter %
% time series as well as the velocity vector timeseries.              %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref = meta_data;
op = zeros(ref(end,1),1);
vel_x = zeros(ref(end,1),1);
vel_y = zeros(ref(end,1),1);
speed = vel_x;
count = zeros(ref(end,1),1);
pos_x = ref(:,4);
pos_y = ref(:,3);
n = 1;  %counter for frame
d = 0;
c = 0; %count of individuals in each frame
%loop for finding centroid at each frame
for i = 1:length(ref(:,1))-1
    cx = ref(i,4);
    cy = ref(i,3);
    vx = ref(i,4) - ref(i,6);
    vy = ref(i,3) - ref(i,5);
    s_i = (sqrt(vx^2 + vy^2));
    if vx == 0 && vy == 0
        i_vx = 0;
        i_vy = 0;
    else
        i_vx = vx / (sqrt(vx^2 + vy^2));
        i_vy = vy / (sqrt(vx^2 + vy^2));
    end
%         if sqrt((cx-960)^2+(cy-544)^2) > 495  %Boundary condition
%             i_vx = nan;
%             i_vy = nan;
%         end
    if ref(i+1,1) == ref(i,1) && ~isnan(vx) && ~isnan(vx)% && sqrt((cx-960)^2+(cy-544)^2) < 410
        vel_x(n) = vel_x(n) + i_vx;
        vel_y(n) = vel_y(n) + i_vy;
        speed(n) = speed(n) + s_i;
        c = c + 1;
    elseif ref(i+1,1) ~= ref(i,1)
        if ~isnan(vx) && ~isnan(vx)
            vel_x(n) = vel_x(n) + i_vx;
            vel_y(n) = vel_y(n) + i_vy;
            speed(n) = speed(n) + s_i;
            c = c + 1;
        end
        vel_x(n) = vel_x(n)/(c);
        vel_y(n) = vel_y(n)/(c);
        speed(n) = speed(n)/c;
        op(n) = sqrt(vel_x(n)^2+vel_y(n)^2);
        count(n) = c;
        n = n + 1;
        c = 0;
    end
end

% op(op == 0) = nan;
% vel_x(vel_x==0) = nan;
% vel_y(vel_y==0) = nan;
vel_x(vel_x>1.5) = nan;
vel_y(vel_y>1.5) = nan;
vel_x(vel_x<-1.5) = nan;
vel_y(vel_y<-1.5) = nan;

figure,
plot(speed(2:length(speed)-1))
% ylim([0 1])

figure,
hist(speed(2:length(speed)-1),100)
% xlim([0 1])

% sum(isnan(op))
%%
% vel = [vel_x, vel_y];
% figure,
% hist3(vel,'CDataMode','auto','FaceColor','interp')
% xlabel('m_x','FontSize',18,'FontWeight','bold')
% ylabel('m_y','FontSize',18,'FontWeight','bold')
% zlabel('Frequency','FontSize',18,'FontWeight','bold')
% colormap(hot)



