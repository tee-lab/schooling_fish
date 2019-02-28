
switch N
    case 15
        N=15;
        Kall=round(linspace(1,N-1,12));
        
        ensemble=200;
        OB=zeros(numel(Kall),ensemble);
        OBm=zeros(numel(Kall),1);
        OBsd=zeros(numel(Kall),1);
        for K1=1:numel(Kall)
            K=Kall(K1);
            disp(K)
            svname=strcat('AveragingGenGAOptimResults_KLdiv_b2_bounds2_N',num2str(N),...
                '_FixdK_',num2str(K),'.mat');
            load(svname)
            
            [fvalmin,indx]=min(scores);
            Xopt1=population(indx,:);
            parfor ll=1:ensemble
                OB(K1,ll)=FishAlignmentModel_AveragingGeneralization(Xopt1,K,data,N,it,0);
            end
        end
        
        for pp=1:numel(Kall)
            [~,temp2]=sort(OB(pp,:));
            OBm(pp)=mean(OB(pp,temp2(1:end/2)));
            OBsd(pp)=std(OB(pp,temp2(1:end/2)));
        end
        
        save('KLdiv_meanSD_N15')
        disp(['KLdiv mean' '   ...   ', 'KL divSD'])
        disp([OBm OBsd])
    case 30
        N=30;
        Kall=round(linspace(1,N-1,12));
        
        ensemble=200;
        OB=zeros(numel(Kall),ensemble);
        OBm=zeros(numel(Kall),1);
        OBsd=zeros(numel(Kall),1);
        for K1=1:numel(Kall)
            K=Kall(K1);
            disp(K)
            svname=strcat('AveragingGenGAOptimResults_KLdiv_b2_bounds2_N',num2str(N),...
                '_FixdK_',num2str(K),'.mat');
            load(svname)
            
            [fvalmin,indx]=min(scores);
            Xopt1=population(indx,:);
            parfor ll=1:ensemble
                OB(K1,ll)=FishAlignmentModel_AveragingGeneralization(Xopt1,K,data,N,it,0);
            end
        end
        
        for pp=1:numel(Kall)
            [~,temp2]=sort(OB(pp,:));
            OBm(pp)=mean(OB(pp,temp2(1:end/2)));
            OBsd(pp)=std(OB(pp,temp2(1:end/2)));
        end
        save('KLdiv_meanSD_N30')
        disp(['KLdiv mean' '   ...   ', 'KL divSD'])
        disp([OBm OBsd])
        
    case 60
        N=60;
        Kall=round(linspace(1,N-1,12));
        
        ensemble=200;
        OB=zeros(numel(Kall),ensemble);
        OBm=zeros(numel(Kall),1);
        OBsd=zeros(numel(Kall),1);
        for K1=1:numel(Kall)
            K=Kall(K1);
            disp(K)
            svname=strcat('AveragingGenGAOptimResults_KLdiv_b2_bounds2_N',num2str(N),...
                '_FixdK_',num2str(K),'.mat');
            load(svname)
            
            [fvalmin,indx]=min(scores);
            Xopt1=population(indx,:);
            parfor ll=1:ensemble
                OB(K1,ll)=FishAlignmentModel_AveragingGeneralization(Xopt1,K,data,N,it,0);
            end
        end
        
        for pp=1:numel(Kall)
            [~,temp2]=sort(OB(pp,:));
            OBm(pp)=mean(OB(pp,temp2(1:end/2)));
            OBsd(pp)=std(OB(pp,temp2(1:end/2)));
        end
        disp(['KLdiv mean' '   ...   ', 'KL divSD'])
        disp([OBm OBsd])
        save('KLdiv_meanSD_N60')
end



