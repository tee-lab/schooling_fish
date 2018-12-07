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


