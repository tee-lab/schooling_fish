%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Jitesh Jhawar                                               %
% Description:                                                        %
%                                                                     %
% This code reads through meta_data and generates the order parameter %
% time series as well as the velocity vectoral components timeseries. %
% Once you run this code, you can either immediately run the dirft    %
% diffusion code that is next in line or save the vel_x.mat and       %
% vel_y.mat to be analyzed later.                                     %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ref = meta_data;
op = zeros(ref(end,1),1);
vel_x = zeros(ref(end,1),1);
vel_y = zeros(ref(end,1),1);
count = zeros(ref(end,1),1);
pos_x = ref(:,4);
pos_y = ref(:,3);
n = 1;  %counter for frame
d = 0;
c = 1; %count of individuals in each frame
%loop for finding centroid at each frame
for i = 1:length(ref(:,1))-1
    cx = ref(i,4);
    cy = ref(i,3);
    vx = ref(i,4) - ref(i,6);
    vy = ref(i,3) - ref(i,5);
    if vx == 0 && vy == 0
        i_vx = 0;
        i_vy = 0;
    else % normalization
        i_vx = vx / (sqrt(vx^2 + vy^2));
        i_vy = vy / (sqrt(vx^2 + vy^2));
    end
    
    if ref(i+1,1) == ref(i,1) && ~isnan(vx) && ~isnan(vx) % && sqrt((cx-960)^2+(cy-544)^2) < 410 % Uncomment last part for boundary effects analysis. Radius of arena = 580 pixels = 90 cms
        vel_x(n) = vel_x(n) + i_vx;
        vel_y(n) = vel_y(n) + i_vy;
        c = c + 1;
    elseif ref(i+1,1) ~= ref(i,1)
        if c >= 10 % threshhold value beyond which we discard data; 5 for n = 15, 10 for n = 30 and 15 for n = 60
            vel_x(n) = vel_x(n)/(c-1); % c-1 because last individual gets skipped in each frame. But that does not affect much.
            vel_y(n) = vel_y(n)/(c-1);
            op(n) = sqrt(vel_x(n)^2+vel_y(n)^2);
        end
        count(n) = c;
        n = n + 1;
        c = 1;
    end
end

op(op == 0) = nan;
vel_x(vel_x==0) = nan;
vel_y(vel_y==0) = nan;
vel_x(vel_x>1.5) = nan;
vel_y(vel_y>1.5) = nan;
vel_x(vel_x<-1.5) = nan;
vel_y(vel_y<-1.5) = nan;

figure,
plot(op)
ylim([0 1])

figure,
hist(op(:,1),100)
% xlim([0 1])

sum(isnan(op))
%% Generate 3D histogram of velocity components
vel = [vel_x, vel_y];
figure,
hist3(vel,'CDataMode','auto','FaceColor','interp')
xlabel('m_x','FontSize',18,'FontWeight','bold')
ylabel('m_y','FontSize',18,'FontWeight','bold')
zlabel('Frequency','FontSize',18,'FontWeight','bold')
colormap(hot)


