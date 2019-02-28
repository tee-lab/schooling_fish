clc
clear

%%
N=15; % SETTING the system size
it=1.2*1e6; % Number of iterations

Gen=150; % Maximum generations
EC=10; % Elite solutions- retained after a generation

name=strcat('hist',num2str(N),'_exp.txt');
data=load(name);

%% PAIRWISE COPYING MODEL
disp('Running Pairwise interaction model')
options = gaoptimset('UseParallel', 'always','EliteCount',EC,...
    'Generations',Gen,'Display','iter');
lbounds=[0,0,0];
ubounds=[1 10 pi];
order=2;
[Xopt2,fval2,~,~,population2,scores2]=...
    ga(@(Xvar) FishAlignmentModel_CopyingGeneralization(Xvar,order,data,N,it,0)...
    ,length(lbounds),[],[],[],[],lbounds,ubounds,[],options);

svname=strcat('PairwiseModel_N_',num2str(N),'_Copying_','.mat');
save(svname)
%% TERN- COPYING MODEL
disp('Running Tern- interaction model')
order=3;
lbounds=zeros(order+1,1);
ubounds=[1 10 2*ones(1,(order-2)) pi];

% CUSTOMIZING INITIAL POPULATION TO GA
initpop3=[rand(50,1) 10*rand(50,1) 2*rand(50,order-2) pi*rand(50,1)];
svname=strcat('PairwiseModel_N_',num2str(N),'_Copying_','.mat'); load(svname,'population2','scores2');
[~,ii]=sort(scores2);
for l=1:25;  initpop3(l,:)=[population2(ii(l),1:end-1) zeros(1,order-2) population2(ii(l),end)]; end
for l=26:40; initpop3(l,:)=[population2(ii(1),1:end-1) zeros(1,order-2) population2(ii(1),end)]; end

options = gaoptimset('UseParallel', 'always','EliteCount',EC,'InitialPopulation',initpop3,'Generations',Gen,'Display','iter');
[Xopt3,fval3,~,~,population3,scores3]=...
    ga(@(Xvar) FishAlignmentModel_CopyingGeneralization(Xvar,order,data,N,it,0)...
    ,order+1,[],[],[],[],lbounds,ubounds,[],options);
svname=strcat('TernModel_N_',num2str(N),'_Copying_','.mat');
save(svname)
%% QUART- COPYING MODEL
disp('Running Quart- interaction model')
order=4;
lbounds=zeros(order+1,1);
ubounds=[1 10 2*ones(1,(order-2)) pi];

% CUSTOMIZING INITIAL POPULATION TO GA
initpop4=[rand(50,1) 10*rand(50,1) 2*rand(50,order-2) pi*rand(50,1)];
svname=strcat('PairwiseModel_N_',num2str(N),'_Copying_','.mat'); load(svname,'population2','scores2');
[~,ii]=sort(scores3);
for l=1:10;  initpop4(l,:)=[population2(ii(l),1:end-1) zeros(1,order-2) population2(ii(l),end)]; end
for l=11:20; initpop4(l,:)=[population2(ii(1),1:end-1) zeros(1,order-2) population2(ii(1),end)]; end
svname=strcat('TernModel_N_',num2str(N),'_Copying_','.mat'); load(svname,'population3','scores3');
[~,ii]=sort(scores2);
for l=21:30;  initpop4(l,:)=[population3(ii(l-19),1:end-1) zeros(1,order-3) population3(ii(l-19),end)]; end
for l=31:40;  initpop4(l,:)=[population3(ii(1),1:end-1) zeros(1,order-3) population3(ii(1),end)]; end

options = gaoptimset('UseParallel', 'always','EliteCount',EC,'InitialPopulation',initpop4,'Generations',Gen,'Display','iter');
[Xopt4,fval4,~,~,population4,scores4]=...
    ga(@(Xvar) FishAlignmentModel_CopyingGeneralization(Xvar,order,data,N,it,0)...
    ,order+1,[],[],[],[],lbounds,ubounds,[],options);
svname=strcat('QuartModel_N_',num2str(N),'_Copying_','.mat');
save(svname)
%% PENTA- COPYING MODEL
disp('Running Penta- interaction model')
order=5;
lbounds=zeros(order+1,1);
ubounds=[1 10 2*ones(1,(order-2)) pi];

% CUSTOMIZING INITIAL POPULATION TO GA
initpop5=[rand(50,1) 10*rand(50,1) 2*rand(50,order-2) pi*rand(50,1)];
svname=strcat('PairwiseModel_N_',num2str(N),'_Copying_','.mat'); load(svname,'population2','scores2');
[~,ii]=sort(scores2);
for l=1:8;  initpop5(l,:)=[population2(ii(l),1:end-1) zeros(1,order-2) population2(ii(l),end)]; end
for l=9:16; initpop5(l,:)=[population2(ii(1),1:end-1) zeros(1,order-2) population2(ii(1),end)]; end
svname=strcat('TernModel_N_',num2str(N),'_Copying_','.mat'); load(svname,'population3','scores3');
[~,ii]=sort(scores3);
for l=17:24;  initpop5(l,:)=[population3(ii(l-15),1:end-1) zeros(1,order-3) population3(ii(l-15),end)]; end
for l=25:32;  initpop5(l,:)=[population3(ii(1),1:end-1) zeros(1,order-3) population3(ii(1),end)]; end
svname=strcat('QuartModel_N_',num2str(N),'_Copying_','.mat'); load(svname,'population4','scores4');
[~,ii]=sort(scores4);
for l=33:40;  initpop5(l,:)=[population4(ii(l-31),1:end-1) zeros(1,order-4) population4(ii(l-31),end)]; end
for l=41:48;  initpop5(l,:)=[population4(ii(1),1:end-1) zeros(1,order-4) population4(ii(1),end)]; end

options = gaoptimset('UseParallel', 'always','EliteCount',EC,'InitialPopulation',initpop5,'Generations',Gen,'Display','iter');
[Xopt5,fval5,~,~,population5,scores5]=...
    ga(@(Xvar) FishAlignmentModel_CopyingGeneralization(Xvar,order,data,N,it,0)...
    ,order+1,[],[],[],[],lbounds,ubounds,[],options);
svname=strcat('PentModel_N_',num2str(N),'_Copying_','.mat');
save(svname)

%% ANALYSIS
disp('AnalysisOfSolutions')
AnalysisOfSolutions