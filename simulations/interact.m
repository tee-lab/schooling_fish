clc;
clear;
close all;
tic
N = 200;
iter = 200000;
x = rand(N,1)*2*pi;
% x = randi([0,1],N,1)*pi;
t = zeros(3*N,1);
pol = zeros(iter,1);
mx = zeros(iter,1);
my = zeros(iter,1);
time = zeros(iter,1);
p = 0.001; s = 0.3; h = 2;% rate parameters
sigma = 3;
Tint = 1/(N*s);
steps = iter;
Tend = steps*Tint;
repetitions = 1;
avgDiff = zeros(102,repetitions);
avgDrift = zeros(102,repetitions);
%assigning clock time to each individual
for j = 1:N
    r1 = rand(1,1);
    t(j) = (1/s)*log(1/r1);
end
for j = (N+1):2*N
    r2 = rand(1,1);
    t(j) = (1/p)*log(1/r2);
end
for j = 2*N+1:3*N
    r3 = rand(1,1);
    t(j) = (1/h)*log(1/r3);
end

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
        %%%%%%% two state case%%%%%%%%
        %             if x(ind) == 0
        %                 x(ind) = pi;
        %             else
        %                 x(ind) = 0;
        %             end
        %%%%%%% two state case%%%%%%%%
        %       x(ind) = ((rand(1,1)*2))*pi;
        t(ind) = (1/s)*log(1/rand(1,1));
        
        % feedback reaction
    elseif ind > N && ind <=2*N
        
        
        t(ind) = (1/p)*log(1/rand(1,1));
        ind = ind - N;
        ind1 = ceil(N*rand(1,1));
        
        %             x(ind) = mode(x);
        x(ind) = x(ind1);   %pairwise interaction
        %             x(ind) = mean(x);
        
        % higher order interaction
    else
        
        t(ind) = (1/h)*log(1/rand(1,1));
        ind = ind - 2*N;
        ind1 = ceil(N*rand(1,1));
        ind2 = ceil(N*rand(1,1));
        %             while ind2 == ind1
        %                 ind2 = ceil(N*rand(1,1));
        %             end
        vec = [sin(x(ind)),cos(x(ind))];
        vec1 = [sin(x(ind1)),cos(x(ind1))];
        vec2 = [sin(x(ind2)),cos(x(ind2))];
        cosTheta1 = dot(vec,vec1)/(norm(vec)*norm(vec1));
        theta1 = acos(cosTheta1);
        cosTheta2 = dot(vec2,vec1)/(norm(vec2)*norm(vec1));
        theta2 = acos(cosTheta2);
        cosTheta3 = dot(vec,vec2)/(norm(vec)*norm(vec2));
        theta3 = acos(cosTheta3);
        if theta1 < theta2 && theta1 < theta3
            c1 = rand(1,1);
            if c1 <= 0.5
                x(ind2) = x(ind);
%                 x(ind1) = x(ind);
            else
                x(ind2) = x(ind1);
%                 x(ind) = x(ind1);
            end
        end
        if theta3 < theta1 && theta3 < theta2
            c1 = rand(1,1);
            if c1 <= 0.5
                x(ind1) = x(ind);
%                 x(ind2) = x(ind);
            else
                x(ind1) = x(ind2);
%                 x(ind) = x(ind2);
            end
        end
        if theta2 < theta1 && theta2 < theta3
            c1 = rand(1,1);
            if c1 <= 0.5
                x(ind) = x(ind1);
%                 x(ind2) = x(ind1);
            else
                x(ind) = x(ind2);
%                 x(ind1) = x(ind2);
            end
        end
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
X = pol.^2;
acf = autocorr(X,2000);
t_lag = 1:2001;
t_lag = t_lag';
f = fit(t_lag(:,1),acf(:,1),'exp1');
val = coeffvalues(f);
b = val(1,2);
ct = abs(ceil(1/(b)));
% [avgDiff(:,1),avgDrift(:,1),op] = driftAndDiffusion(X,Tint,40,N,p,s,sigma,h);
figure,
histogram(pol,40,'Normalization','probability')
title(['N = ' num2str(N) ' f = ' num2str(p) ' h = ' num2str(h) ' s = ' num2str(s) ' sigma = pi/' num2str(sigma)])
xlabel('Group alignment','fontsize',14,'FontWeight','bold')
ylabel('Probability','fontsize',14,'FontWeight','bold')
xlim([0,1])

vel_y = my;
vel_x = mx;
save('interact_vel_x_y.mat','vel_x', 'vel_y');
toc
