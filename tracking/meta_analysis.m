%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Amith Kumar U R                                             %
% Description:                                                        %
%                                                                     %
% This code stitches all the tracked data generated using             %
% loop_track.m. For input just select the folder in which all         %
% the data is saved by loop_track.m. The output is a file called      %
% trac_metadata.mat that can be used by op_calculate.m to calculate   %
% polarization.                                                       %
%                                                                     %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; 
%%
mdir =uigetdir;
fils_path = [mdir '/trac*'];
files=dir(fils_path) ; 
%%
for i=1:size(files,1) 
    
    a{i} = str2mat(files(i).name);
end
b = sort_nat(a);
%%
for i=1:length(b)
%     clear sept;
    dtid{1,i}=load(b{i});
%     Imgs{1,i}=imread(files(i).name);
end
%%
clear meta_data;
d = dtid{1,1}.data_sheet;
for k=2:size(dtid,2)
    clear d1 d2 d2m 
   
    d1 = dtid{1,k}.data_sheet;
    d2 = d1(~isnan(d1(:,1)),:);
    d2 = d2(d2(:,1) > 31,:);
    d2m = d2(:,1) - 31 + max(d(:,1));
    d2(:,1) = d2m;
    d = cat(1,d,d2);
end
meta_data = d;
clear d;
%%
count = zeros(max(meta_data(:,1)),3);
parfor l=2:max(meta_data(:,1))
    test_frame = meta_data(meta_data(:,1) == l,:);
    count1(l,1) = size(test_frame,1);
    count2(l,1) = size(test_frame(~isnan(test_frame(:,5)),:),1);
    testf1 = sqrt(test_frame(:,7).^2 + test_frame(:,8).^2);
    count3(l,1) = size(testf1(testf1 ~= 0,:),1);
    l
end
count = cat(2,count1,count2,count3);
%%
clear filtop_p filtop_r ;
p=10;
% [op,oplocal] = New_OP_polarization(meta_data,1,max(meta_data(:,1)),100,count,30);
%[filtop_op,filtop_local] = New_OP_polarization(meta_data,1,max(meta_data(:,1)),30,count,p);
[filtop_p,filtop_r] = New_OP(meta_data,1,max(meta_data(:,1)),count,p);
t=(1:max(meta_data(:,1)))/(25*60);

save([fils_path(1:end-1) '_metadata.mat'],'meta_data','filtop_p','filtop_r','count','t')
% mindistance;
% %% Cluster analysis
% clear cd ctoff;
% cd = cell(max(meta_data(:,1)),1);
% ctoff = zeros(max(meta_data(:,1)),1);
% 
% parfor i=2:max(meta_data(:,1))
%     
%     ptest= meta_data(meta_data(:,1) == i, :);
%         if size(ptest,1) > 20
%             if size(ptest,1) < 5
%                 kmax = size(ptest,1);
%             else
%                 kmax = 5;
%             end
% 
%         eva = evalclusters(ptest(:,3:4),'kmeans','gap','klist',[1:kmax],'SearchMethod','firstMaxSE');
%         ctoff(i,1) = eva.OptimalK;
%     
%             if ctoff(i,1) > 1
%                 cd{i,1} = clusterdata(ptest(:,3:4),'maxclust',ctoff(i,1));
%             else
%                 cd{i,1} = ones(size(ptest,1),1);
%             end
%         else
%         cd{i,1} = nan;
%         ctoff(i,1) = nan;
%     end
%     i
% end
% 
% clear grpid
% grpid=cell(length(cd),5);
% 
% for i=1:length(cd)
%     if ~isnan(cd{i,1})
%         dm = cd{i,1};
%         ptest= meta_data(meta_data(:,1) == i, :);
%         mx = max(dm);
%         for j=1:mx
%             grpid{i,j} = ptest(dm==j,:);
%         end
%     else
%         grpid{i} =nan;
%     end
%     
% end
% 
% fulname = [fils_path(1:end-1) 'cluster_data.mat'];
% % path = mfile.Path;
% % Name = videor.Name(1:end-4) ; 
% save(fulname, 'grpid','cd','ctoff');
% %%
% % w=24; %Window size for velocity. 600 = 1 minute.
% % [meanvel,meanve]=velocitis(dtid,w);
% % %%
% % %tic;[distance_covers]=distancs(dtid);toc
% % %%
% % tic;
% % for k=1:length(dtid)
% %       clear m n
% %       [m]=nearestnghbr_dist(dtid{1,k}.data_sheet);  % use "nearestnghbr_distance" for range wise nbr dist
% %       nbrdistance{1,k}=m;   % mean neighbor distance in each frame
% %       nbrdistance{2,k}=mean(m);
% %     %   frame_avg_nbrdist{1,k}=n;% neighbor distance 
% %     %   frame_avg_nbrdist{2,k}=mean(n);
% %     %   
% %       k
% % end
% % % toc;
% % 
% % clear nbrdist;
% % l=1;
% % for i=1:length(nbrdistance)
% %     w=24;
% %     for x=1:w:length(nbrdistance{1,i})-w
% %         nbrdist(1,l)=mean(nbrdistance{1,i}(x:x+(w-1),1));
% %         % nbrdist(2,l)=mean(nbrdistance{1,i}(x:x+(w-1),2));
% %         % nbrdist(3,l)=mean(nbrdistance{1,i}(x:x+(w-1),3));
% %         nbrdist(2,l)=std(nbrdistance{1,i}(x:x+(w-1),1));
% %         nbrdist(3,l)=var(nbrdistance{1,i}(x:x+(w-1),1));
% %         l=l+1;
% %     end
% % end
% % toc
