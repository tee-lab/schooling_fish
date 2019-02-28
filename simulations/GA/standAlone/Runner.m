clc
clear

%%
N=15; % System size
iter=1e5; % Timesteps
kk=2;

name=strcat('hist',num2str(N),'_exp.txt');
data=load(name);

% For averaging interactions
alneigh=2;
Xvar=[0.1 2 pi]; % Size is 3x1: [spontaneous, pairwise, fluc]
[Ob1,polx1, poly1, pol1, t1]=FishAlignmentModel_AveragingGeneralization(Xvar,alneigh,data,N,iter,kk);

% For copying interactions
order=3;
Xvar=[0.1 2 0.1 pi]; % Size should be (order+1): [spontaneous, pairwise, ternary, fluc]
[Ob2,polx2, poly2, pol2, t2]=FishAlignmentModel_CopyingGeneralization(Xvar,order,data,N,iter,kk);