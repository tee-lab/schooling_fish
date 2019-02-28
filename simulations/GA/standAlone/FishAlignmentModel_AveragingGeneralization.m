function [Ob,polx, poly, pol, t]=JithVishmodel_sVicsek_fixdK(Xvar,alneigh,data,N,iter,kk)


spon=Xvar(1);tur=Xvar(3);
Vic=Xvar(2);

% N = 15;
% iter = 5*1e5;
% dt = 1/(N*spon);

xs = rand(N,1)*2*pi; % State of the system
cl_time = zeros(2*N,1); % Clock time

pol = zeros(iter,1); % Order parameter- polarization
polx = zeros(iter,1); % Order parameter- polarization
poly = zeros(iter,1); % Order parameter- polarization
t = zeros(iter,1); % Event times

% Tend = iter*dt;
% repetitions = 1; % Independent realizations

% Assigning clock time to each individual
cl_time(1:N)=1./(spon).*log(1./(rand(N,1)));
cl_time(N+1:2*N)=1./(Vic).*log(1./(rand(N,1)));

%% TIME EVOLUTION OF STATES
for i=1:iter
    ind = find(cl_time==min(cl_time),1);
    t(i) = cl_time(ind);
    
    %     Updated clock
    cl_time = cl_time - cl_time(ind);
    
    if ind <= N %Spontaneous change
        
        ra = normrnd(0,tur);  %Gaussian noise
        while (ra < -pi || ra > pi)
            ra = normrnd(0,tur);
        end
        
        xs(ind) = rem(xs(ind) + ra,2*pi); % Getting it between 0 and 2*pi
        if xs(ind)<0
            xs(ind)=xs(ind)+2*pi;
        end
        
        cl_time(ind) = (1/spon)*log(1/rand(1,1));
        
    else %Vicsek alignement interaction
        
        cl_time(ind) = (1/Vic)*log(1/rand(1,1));
        ind = ind - N;
        
        neighall=1:N;
        neighall(ind)=[];
        neighs=neighall(randperm(N-1,alneigh));
        
        xs(ind)=atan2(mean(sin(xs(neighs)),1),mean(cos(xs(neighs)),1));
        if xs(ind)<0
            xs(ind)=2*pi+xs(ind);
        end
    end
    vxm = mean(cos(xs));
    vym = mean(sin(xs));
    polx(i,1) = vxm;
    poly(i,1) = vym;
    pol(i,1) = norm([vxm;vym]);
end
%% Calculating objective function value
[freq]=hist(pol,data(:,1));
P=data(:,2)./sum(data(:,2));
Q=(freq./sum(freq))';

% LeastSquares
% Ob=sqrt(sum((P-Q).^2));

% KLDIV
KL=Q.*log2(Q./P);
KL(Q==0)=0;
Ob=sum(KL);

switch kk
    case 0
    case 1
        figure
        plot(data(:,1),Q)
        hold all
    case 2
        figure
        bar(data(:,1),P); hold all; alpha(0.7)
        histogram(pol,40,'Normalization','probability')
        alpha(0.4)
        legend('data','model')
        hold all
end
end
