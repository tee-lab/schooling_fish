clc
clear
%%
N=15; % System size
it=5e5; % Timesteps
Gen=150;
name=strcat('hist',num2str(N),'_exp.txt');
data=load(name);

%% GA simulation
lbounds=[0,0,0]; ubounds=[10,10,pi];
options = gaoptimset('UseParallel', 'always','Display','iter','Generations',Gen);

for K=round(linspace(1,N-1,12)) % Runs for 15 cases, equally distributed between 1 and N
    
    [Xopt,fval,exitflag,output,population,scores]=...
        ga(@(Xvar) FishAlignmentModel_AveragingGeneralization(Xvar,K,data,N,it,0)...
        ,length(lbounds),[],[],[],[],lbounds,ubounds,[],options);
    
    svname=strcat('AveragingGenGAOptimResults_KLdiv_b2_bounds2_N',num2str(N),'_FixdK_',num2str(K),'.mat');
    save(svname)s
end

%% Analysis of the optimized solutions
AnalysisOfOptimizedsolns