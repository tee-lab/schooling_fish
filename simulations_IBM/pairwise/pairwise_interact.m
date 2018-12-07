% clear;
N = 30;
iter = 1000;
x = rand(N,1)*2*pi;
t = zeros(2*N,1);
pol = zeros(iter,1);
mx = zeros(iter,1);
my = zeros(iter,1);
time = zeros(iter,1);
r = 3.7; s = 0.015; % rate parameters
eta = 1;
sigma = 3;
Tint = 0.12;
steps = iter;
Tend = steps*Tint;
repetitions = 1;
avgDiff = zeros(102,repetitions);
avgDrift = zeros(102,repetitions);
% for rep = 1:repetitions
%assigning clock time to each individual
for j = 1:N
    r1 = rand(1,1);
    t(j) = (1/s)*log(1/r1);
end
for j = (N+1):2*N
    r2 = rand(1,1);
    t(j) = (1/r)*log(1/r2);
end
c = 0;


%begin iteration
T = 0;
Tprint = 0;
i = 1;
n = 1;
while (T < Tend)
    
    %find minimum clock time
    ind = find(t==min(t));
    time(i) = t(ind);
    T = T + t(ind);
    t = t - t(ind);
    % spontaneous reaction
    if ind <= N
        
        ra = normrnd(0,pi/sigma);  %Gaussian noise
        while (ra < -pi || ra > pi)
            ra = normrnd(0,pi/sigma);
        end
        
        x(ind) = x(ind) + ra;   %Gaussian noise
        
        %       x(ind) = ((rand(1,1)*2))*pi;
        t(ind) = (1/s)*log(1/rand(1,1));
        
    else
        % feedback reaction
        t(ind) = (1/r)*log(1/rand(1,1));
        ind = ind - N;
        ind1 = ceil(N*rand(1,1));
        
        %         x(ind) = mode(x);
        x(ind) = x(ind1);   %pairwise interaction
        %                 x(ind) = mean(x);
    end
    vx = cos(x);
    vy = sin(x);
    
    if T > Tprint
        pol(n,1) = (((sum(vx))^2+(sum(vy))^2)^0.5)/N;
        mx(n,1) = mean(vx);
        my(n,1) = mean(vy);
        n = n + 1;
        Tprint = Tprint + Tint;
        
    end
    
    i = i + 1;
end

% figure,
% plot(pol)
% ylim([0,1])
% figure,
% hist(mx,100)
% xlim([-1,1])
% title(['N = ' num2str(N) ' f = ' num2str(r) ' s = ' num2str(s)])
X = pol.^2;

% figure,
% hist(my,100)
% xlim([-1,1])

%%
acf = autocorr(mx,250);
t_lag = 1:251;
t_lag = t_lag';
f = fit(t_lag(:,1),acf(:,1),'exp1');
val = coeffvalues(f);
b = val(1,2);
ct = abs(ceil(1/(b)));
%
% figure,
% plot(t_lag,acf)
% xlim([0 250])
% xlabel('Time Lag')
% ylabel('ACF')
% title(['N = ' num2str(N) ' f = ' num2str(r) ' s = ' num2str(s)])
% hline = refline([0 0]);
% hline.Color = 'r';



% end
figure,
histogram(pol,40,'Normalization','probability')
title(['N = ' num2str(N) ' f = ' num2str(r) ' s = ' num2str(s) ' sigma = pi/' num2str(sigma)])
xlabel('Group alignment','fontsize',14,'FontWeight','bold')
ylabel('Probability','fontsize',14,'FontWeight','bold')
xlim([0,1])
%%
% [avgDiff(:,1),avgDrift(:,1),op] = driftAndDiffusion(X,Tint,40,N,r,s,sigma);
vel_y = my;
vel_x = mx;
