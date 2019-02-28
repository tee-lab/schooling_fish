% AIC analysis
ensemble=200;
kk=0;samples=40;

OB=zeros(4,ensemble,1);
OBm=zeros(4,1);
AIC=zeros(4,ensemble,1);
RSS=zeros(4,ensemble,1);
AICm=zeros(4,1);
AICsd=zeros(4,1);
OBsd=zeros(4,1);

load('PairwiseModel_N_60_Copying_.mat')
[~,ii]=sort(scores2);
Xvar=population2(ii(1),:);
order=2;
parfor pp=1:ensemble
    [ob,Rss]= FishAlignmentModel_CopyingGeneralizationRSS(Xvar,order,data,N,it,kk);
    RSS(1,pp,1)=Rss;
    AIC(1,pp,1)=samples*log(Rss/samples)+2*(order+1);
    OB(1,pp,1)=ob;
end

load('TernModel_N_60_Copying_.mat')
[~,ii]=sort(scores3);
Xvar=population3(ii(1),:);
parfor pp=1:ensemble
    [ob,Rss]= FishAlignmentModel_CopyingGeneralizationRSS(Xvar,order,data,N,it,kk);
    RSS(2,pp,1)=Rss;
    AIC(2,pp,1)=samples*log(Rss/samples)+2*(order+1);
    OB(2,pp,1)=ob;
end

load('QuartModel_N_60_Copying_.mat')
[~,ii]=sort(scores4);
Xvar=population4(ii(1),:);
parfor pp=1:ensemble
    [ob,Rss]= FishAlignmentModel_CopyingGeneralizationRSS(Xvar,order,data,N,it,kk);
    RSS(3,pp,1)=Rss;
    AIC(3,pp,1)=samples*log(Rss/samples)+2*(order+1);
    OB(3,pp,1)=ob;
end

load('PentModel_N_60_Copying_.mat')
[~,ii]=sort(scores5);
Xvar=population5(ii(1),:);
parfor pp=1:ensemble
    [ob,Rss]= FishAlignmentModel_CopyingGeneralizationRSS(Xvar,order,data,N,it,kk);
    RSS(4,pp,1)=Rss;
    AIC(4,pp,1)=samples*log(Rss/samples)+2*(order+1);
    OB(4,pp,1)=ob;
end

%% IDENTIFYING CANDIDATE RUNS
for pp=1:order-1
    [~,temp2]=sort(OB(pp,:,1));
    AICm(pp,1)=mean(AIC(pp,temp2(1:end/2)));
    AICsd(pp,1)=std(AIC(pp,temp2(1:end/2)));
    OBm(pp,1)=mean(OB(pp,temp2(1:end/2)));
    OBsd(pp,1)=std(OB(pp,temp2(1:end/2)));
end
disp(['OBm' ' ...  ' 'OBsd' '  ... ' 'AICm' ' ...  ' 'AICsd'])
disp([OBm OBsd AICm AICsd])
