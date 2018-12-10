%% create sequence of vectoral order parameter
inc_x = 0.1;
inc_y = 0.1;
op_x = -1:inc_x:1;
op_y = -1:inc_y:1;
%%
n = size(vel_y); %store number of data points
diffY = zeros(n); %Diffusion function of vel_y
diffXY = zeros(n);  %Cross diffusion function of vel_x and vel_y
driftY = zeros(n);  %Drift function of vel_y
avgDriY = zeros(size(op_y,2),size(op_x,2));  %Stores the binned Drift Y function 
avgDifY = zeros(size(op_y,2),size(op_x,2));  %Stores the binned Diffusion Y function
avgDifXY = zeros(size(op_y,2),size(op_x,2)); %Stores the binned cross diffusion function

diffX = zeros(n);  %Diffusion function of vel_x
driftX = zeros(n);  %Drift function of vel_x
avgDriX = zeros(size(op_x,2),size(op_y,2));  %Stores the binned Drift X function 
avgDifX = zeros(size(op_x,2),size(op_y,2));  %Stores the binned Difusion X function

%% lets calculate drift and diffusion for each data point
dt = 0.12;  %change this to 1/(Ns) if analysing data from ternary interaction model. N is system size and s is spontaneous reaction rate 
Dt = 10;
for i = 1:size(vel_y) - Dt
    diffY(i) = ((vel_y(i+1) - vel_y(i))^2) / dt;
    diffXY(i) = ((vel_y(i+1) - vel_y(i))*(vel_x(i+1) - vel_x(i))) / dt;
    driftY(i) = (vel_y(i+Dt) - vel_y(i)) / (Dt*dt);
end
for i = 1:(size(vel_x) - Dt)
    diffX(i) = ((vel_x(i+1) - vel_x(i))^2) / dt;
    driftX(i) = (vel_x(i+Dt) - vel_x(i)) / (Dt*dt);
end

%% lets start binning the data for Drift Y and Diffusion Y function
bin_x = -1;
m = 1;

while (bin_x < 1)
    n = 1;
    bin_y = -1;
    while (bin_y < 1)
        c1 = 0;
        c2 = 0;
        binDriY = zeros(length(vel_y),1);
        binDifY = zeros(length(vel_y),1);
        for i = 1:size(vel_y)
            if (vel_x(i) < bin_x + inc_x && vel_x(i) >= bin_x && vel_y(i) < bin_y + inc_y && vel_y(i) >= bin_y && ~isnan(driftY(i)))
                avgDriY(n,m) = avgDriY(n,m) + driftY(i);
                c1 = c1 + 1;
            end
            if (vel_x(i) < bin_x + inc_x && vel_x(i) >= bin_x && vel_y(i) < bin_y + inc_y && vel_y(i) >= bin_y && ~isnan(diffY(i)))
                avgDifY(n,m) = avgDifY(n,m) + diffY(i);
                avgDifXY(n,m) = avgDifXY(n,m) + diffXY(i);
                c2 = c2 + 1;
            end
            %             binDriY(j) = driftY(i);
            %             binDifY(j) = diffY(i);
            %             j = j + 1;
        end
        
        %         binDriY(binDriY==0) = nan;
        %         binDifY(binDifY==0) = nan;
        if c1 > 0
            avgDriY(n,m) = avgDriY(n,m)/c1;
        end
        if c2 > 0
            avgDifY(n,m) = avgDifY(n,m)/c2;
            avgDifXY(n,m) = avgDifXY(n,m)/c2;
        end
        
        bin_y = bin_y + inc_y;
        n = n + 1;
    end
    bin_x = bin_x + inc_x;
    m = m + 1;
end
%% binning the data for Drift X and Diffusion X function
bin_x = -1;
m = 1;

while (bin_x < 1)
    n = 1;
    bin_y = -1;
    while (bin_y < 1)
        c1 = 0;
        c2 = 0;
        binDriX = zeros(length(vel_x),1);
        binDifX = zeros(length(vel_x),1);
        for i = 1:size(vel_x)
            if (vel_x(i) < bin_x + inc_x && vel_x(i) >= bin_x && vel_y(i) < bin_y + inc_y && vel_y(i) >= bin_y && ~isnan(driftY(i)))
                avgDriX(n,m) = avgDriX(n,m) + driftX(i);
                c1 = c1 + 1;
            end
            if (vel_x(i) < bin_x + inc_x && vel_x(i) >= bin_x && vel_y(i) < bin_y + inc_y && vel_y(i) >= bin_y && ~isnan(diffY(i)))
                avgDifX(n,m) = avgDifX(n,m) + diffX(i);
                c2 = c2 + 1;
            end
            %             binDriY(j) = driftY(i);
            %             binDifY(j) = diffY(i);
            %             j = j + 1;
        end
        
        %         binDriY(binDriY==0) = nan;
        %         binDifY(binDifY==0) = nan;
        if c1 > 0
            avgDriX(n,m) = avgDriX(n,m)/c1;
        end
        if c2 > 0
            avgDifX(n,m) = avgDifX(n,m)/c2;
        end
        
        bin_y = bin_y + inc_y;
        n = n + 1;
    end
    bin_x = bin_x + inc_x;
    m = m + 1;
end

%% Plotting begins
%Time Series
% figure
% subplot(2,1,1)       % add first plot in 2 x 1 grid
% scatter(time(1:4000),vel_x(1:4000))
% xlabel('Frame','FontSize',15,'FontWeight','bold')
% ylabel('Polarization (X)','FontSize',12,'FontWeight','bold')
% ylim([-1 1])
% % histogram
% subplot(2,1,2)       % add second plot in 2 x 1 grid
% histogram(vel_x,50)
% ylabel('Frequency','FontSize',12,'FontWeight','bold')
% xlabel('Polarization (X)','FontSize',15,'FontWeight','bold')
% xlim([-1, 1])
%% ploting drift Y with polarization X
%subplot(4,1,3)
% ind = 14
% figure,
% scatter(op_x,avgDriY(ind,:))
% title(['Vy = ' num2str(op_y(ind))])
% ylabel('Deterministic Factor (Y)','FontSize',15,'FontWeight','bold')
% xlabel('Polarization (X)','FontSize',15,'FontWeight','bold')
% ylim([-0.006, 0.006])
%% ploting diffusion Y with polarization X
%subplot(4,1,4)
% ind = 18
% figure,
% scatter(op_x,avgDifY(ind,:))
% title(['Vy = ' num2str(op_y(ind))])
% xlabel('Polarization (X)','FontSize',15,'FontWeight','bold')
% ylabel('Stochastic Factor (Y)','FontSize',15,'FontWeight','bold')
% ylim([0, 0.002])
%% plotting drift Y with polarization Y
% ind = 17
% figure,
% scatter(op_y,avgDriY(:,ind))
% title(['Vx = ' num2str(op_x(ind))])
% xlabel('Polarization (Y)','FontSize',15,'FontWeight','bold')
% ylabel('Deterministic Factor (Y)','FontSize',12,'FontWeight','bold')
% ylim([-0.005, 0.005])
% %% ploting diffusion Y with polarization Y
% ind = 17
% figure,
% scatter(op_y,avgDifY(:,ind))
% title(['Vx = ' num2str(op_x(ind))])
% xlabel('Polarization (Y)','FontSize',15,'FontWeight','bold')
% ylabel('Stochastic Factor (Y)','FontSize',12,'FontWeight','bold')
% ylim([0, 0.004])
%% 3D plotting Drift
% figure,
% surf(op_x,op_y,avgDriY)             % The surface plotting function.
% xlabel('Polarization (X)','FontSize',15,'FontWeight','bold')
% ylabel('Polarization (Y)','FontSize',15,'FontWeight','bold')
% zlabel('Deterministic Factor (Y)','FontSize',15,'FontWeight','bold')
% h = colorbar
% h.Limits = [-0.06 0.03]
%% 3D plotting Diffusion
% figure,
% surf(op_x,op_y,avgDifY)             % The surface plotting function.
% xlabel('Polarization (X)','FontSize',15,'FontWeight','bold')
% ylabel('Polarization (Y)','FontSize',15,'FontWeight','bold')
% zlabel('Stochastic Factor (Y)','FontSize',15,'FontWeight','bold')
% h = colorbar
% h.Limits = [0 0.005]
%%
vel = [vel_x, vel_y];
figure,
hist3(vel,'CDataMode','auto','FaceColor','interp')
xlabel('m_x','FontSize',18,'FontWeight','bold')
ylabel('m_y','FontSize',18,'FontWeight','bold')
zlabel('Frequency','FontSize',18,'FontWeight','bold')
colormap(hot)
%%
figure,
a = repmat(op_x,length(op_x),1);
b = repmat(op_y,length(op_y),1);
a = a(:); b = b'; b = b(:);
avgDifY(avgDifY==0) = nan;
scatter3(a(:),b(:),avgDifY(:),'filled','red')             % The surface plotting function.
% hold on
% scatter3(a(:),b(:),avgDifY30(:))
% hold on
% scatter3(a(:),b(:),avgDifY60(:))
xlabel('m_x','FontSize',15,'FontWeight','bold')
ylabel('m_y','FontSize',15,'FontWeight','bold')
zlabel('Stochastic Factor g^2_m_y(m_x,m_y)','FontSize',15,'FontWeight','bold')
legend('N = 15', 'N = 30','N = 60','Location','north')
%%
figure,
a = repmat(op_x,length(op_x),1);
b = repmat(op_y,length(op_y),1);
a = a(:); b = b'; b = b(:);
avgDifX(avgDifX==0) = nan;
scatter3(a(:),b(:),avgDifX(:),'filled','red')             % The surface plotting function.
% hold on
% scatter3(a(:),b(:),avgDifY30(:))
% hold on
% scatter3(a(:),b(:),avgDifY60(:))
xlabel('m_x','FontSize',15,'FontWeight','bold')
ylabel('m_y','FontSize',15,'FontWeight','bold')
zlabel('Stochastic Factor g^2_m_x(m_x,m_y)','FontSize',15,'FontWeight','bold')
legend('N = 15', 'N = 30','N = 60','Location','north')
% zlim([0,0.5])

%% Overlaying different Group size results. You may save drift and diffusion function for 30 and 60 first as their name for eg. avgDriX30.
figure,
a = repmat(op_x,length(op_x),1);
b = repmat(op_y,length(op_y),1);
a = a(:); b = b'; b = b(:);
S = repmat(ones(length(op_x),1)*30,numel(op_x),1);
s = S(:);
C = repmat((1:5:5*length(op_x)),numel(op_x),1);
c = C(:);
avgDriY(avgDriY==0) = nan;
scatter3(a(:),b(:),avgDriY(:),'filled','red');%s,c);             % The surface plotting function.
% hold on
% scatter3(a(:),b(:),avgDriY30(:),'filled');%s,c);
% hold on
% scatter3(a(:),b(:),avgDriY60(:),'filled');%s,c);
xlabel('m_x','FontSize',15,'FontWeight','bold')
ylabel('m_y','FontSize',15,'FontWeight','bold')
zlabel('Deterministic Factor f_m_y(m_x,m_y)','FontSize',15,'FontWeight','bold')
legend('N = 15', 'N = 30','N = 60','Location','north')

%%
figure,
a = repmat(op_x,length(op_x),1);
b = repmat(op_y,length(op_y),1);
a = a(:); b = b'; b = b(:);
S = repmat(ones(length(op_x),1)*30,numel(op_x),1);
s = S(:);
C = repmat((1:5:5*length(op_x)),numel(op_x),1);
c = C(:);
avgDriX(avgDriX==0) = nan;
scatter3(a(:),b(:),avgDriX(:),'filled','red');%s,c);             % The surface plotting function.
% hold on
% scatter3(a(:),b(:),avgDriX30(:),'filled');%s,c);
% hold on
% scatter3(a(:),b(:),avgDriX60(:),'filled');%s,c);
xlabel('m_x','FontSize',15,'FontWeight','bold')
ylabel('m_y','FontSize',15,'FontWeight','bold')
zlabel('Deterministic Factor f_m_x(m_x,m_y)','FontSize',15,'FontWeight','bold')
legend('N = 15', 'N = 30','N = 60','Location','north')
%% fitting
% f = fit([a,b],avgDriY(:),'poly23')
% plot(f,[a,b],avgDriY(:))

%% Fitting of various functional forms to Drift and Diffusion can be Done here
xyz = [a b avgDriY(:)]; %x y z are column vextors; you can change avgDriY to avgDriX or avgDifX or avgDifY
xyz_no_nans = xyz(~any(isnan(xyz),2),:); %keep rows which have no nans
x = xyz_no_nans(:,1); %etc
y = xyz_no_nans(:,2);
z = xyz_no_nans(:,3);
X0 = x;
X1 = [x,y];
X2 = [x.^2,y.^2];
X = [x,y,x.*y,x.^2,y.^2,x.^3,y.^3];
md2 = fitlm(X,z)
%% Plotting the fitted function
% figure,
% plot(md1,[x,y],z)
%%
% [f,gof] = fit([x,y],z,'poly23')
% plot(f,[x,y],z)
% xlabel('Polarization (X)','FontSize',15,'FontWeight','bold')
% ylabel('Polarization (Y)','FontSize',15,'FontWeight','bold')
% zlabel('Deterministic Factor f(Y)','FontSize',15,'FontWeight','bold')

%% Overlaying results again
% figure,
% a = repmat(op_x,length(op_x),1);
% b = repmat(op_y,length(op_y),1);
% a = a(:); b = b'; b = b(:);
% avgDifXY(avgDifXY==0) = nan;
% scatter3(a(:),b(:),avgDifXY(:))             % The surface plotting function.
% % hold on
% % scatter3(a(:),b(:),avgDifY30(:))
% % hold on
% % scatter3(a(:),b(:),avgDifY60(:))
% xlabel('m_x','FontSize',15,'FontWeight','bold')
% ylabel('m_y','FontSize',15,'FontWeight','bold')
% zlabel('Stochastic Factor g^2_m_x_y(m_x,m_y)','FontSize',15,'FontWeight','bold')
% legend('N = 30', 'N = 30','N = 60','Location','north')
% zlim([0,0.5])
%% plot time series
% t = (1:length(vel_x))*0.12;
% figure,
% scatter(t,vel_x,'red')
% xlim([4200 5000])
% ylim([-1 1])
% xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
% ylabel('m_x','FontSize',15,'FontWeight','bold')
% 
%% You can plot single vector time series between given time points here
% figure,
% scatter(t,vel_y,'red')
% xlim([4200 4500])
% ylim([-1 1])
% xlabel('Time (seconds)','FontSize',15,'FontWeight','bold')
% ylabel('m_y','FontSize',15,'FontWeight','bold')
