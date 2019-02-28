function [Ob,polx, poly, pol, t]=JithVishmodel_generalCumulativeHigherOrder(Xvar,alneigh,data,N,iter,kk)


spon=Xvar(1);
tur=Xvar(end);

HoI(:,1)=Xvar(2:end-1);

% N = 15;
% iter = 5*1e5;
% dt = 1/(N*spon);

xs = rand(N,1)*2*pi; % State of the system

pol = zeros(iter,1); % Order parameter- polarization
polx = zeros(iter,1); % Order parameter- polarization
poly = zeros(iter,1); % Order parameter- polarization
t = zeros(iter,1); % Event times

% Tend = iter*dt;
% repetitions = 1; % Independent realizations

% Assigning clock time to each individual
cl_time(1)=1./(spon).*log(1./(rand(1,1)));

cl_time(2:alneigh)=1./(HoI).*log(1./(rand(alneigh-1,1)));

%% TIME EVOLUTION OF STATES
for i=1:iter
    ind = find(cl_time==min(cl_time));
    t(i) = cl_time(ind);
    
    %     Updated clock
    cl_time = cl_time - cl_time(ind);
    
    if ind==1
        cl_time(ind) = (1/spon)*log(1/rand(1,1)); %New clock time
        
        ra = normrnd(0,tur);  %Gaussian noise
        while (ra < -pi || ra > pi)
            ra = normrnd(0,tur);
        end
        
        foc=randperm(N,1);
        xs(foc) = rem(xs(foc) + ra,2*pi);
        % Getting it between 0 and 2*pi
        if xs(foc)<0
            xs(foc)=xs(foc)+2*pi;
        end
        
    else
        
        for ev=2:alneigh
            if ind==ev
                cl_time(ind) = (1/HoI(ev-1))*log(1/rand(1,1));
                
                neighs_ind=randperm(N,ev); % Randomly choose a alneigh neighbors
                % Finding the odd one out
                v_neigh=zeros(ev,1);
                for j=1:ev
                    v_neigh(j)=norm([mean(cos(xs(neighs_ind(neighs_ind ~=neighs_ind(j)))))...
                        mean(sin(xs(neighs_ind(neighs_ind ~=neighs_ind(j)))))]);
                end
                [~,foc_ind]=max(v_neigh);
                foc=neighs_ind(foc_ind); % ODD one out
                neighs=neighs_ind(1:end ~=foc_ind); % subset of more aligned neighbors
                
                xs(foc)=xs(neighs(randperm(ev-1,1)));
                break;
            end
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