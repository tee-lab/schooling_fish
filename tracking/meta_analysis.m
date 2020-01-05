%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Amith Kumar U R                                              %
% Description:                                                         %
%                                                                      %
% This code stitches all the tracked data generated using              %
% loop_track.m. Upon running the code you will be asked to select the  %
% folder in which all the data is saved. The output is a file called   %
% trac_rawdata.mat that can be used as an input by calc_polarization.m %
% to calculate polariztion.                                            %
%                                                                      %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc;
clear;
close all;
%%
mdir =uigetdir;
addpath(mdir);
files_path = [mdir '/trac*'];
files=dir(files_path) ; 
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
clear raw_data;
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
raw_data = d;
clear d;
%%
count = zeros(max(raw_data(:,1)),3);
parfor l=2:max(raw_data(:,1))
    test_frame = raw_data(raw_data(:,1) == l,:);
    count1(l,1) = size(test_frame,1);
    count2(l,1) = size(test_frame(~isnan(test_frame(:,5)),:),1);
    testf1 = sqrt(test_frame(:,7).^2 + test_frame(:,8).^2);
    count3(l,1) = size(testf1(testf1 ~= 0,:),1);
    l
end
count = cat(2,count1,count2,count3);

save([mdir '/rawdata_track.mat'],'raw_data');
csvwrite([mdir '/rawdata_track.csv'],raw_data);