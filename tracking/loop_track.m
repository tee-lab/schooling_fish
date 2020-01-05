%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Amith Kumar U R                                                %
% Description:                                                           %
%                                                                        %
% The input for this code is the output file (in .mat format) generated  %
% by detections.m, which has to be selected through the pop-up window.   %
% This code calls the fish_tracker code developed by Student Dave and    %
% can be found at his site http://studentdavestutorials.weebly.com/.     %
% The tracker code uses uses the Kalman-filter method for tracking       %
% individuals between frames. This method generates huge amount of       %
% data as the algorithm progresses and therefore storing it in one       %
% single file makes the running very slow. Therefore, we divide the      %
% data into several chunks by calling fish_tracker multiple times        %
% using this code. Once fish_tracker works on each chunk and we create   %
% a new file to store data at every iteration. Later we stitch all the   %
% data together using meta_analysis.m. The output files calles           %
% tracked***.mat will be saved in the same folder where the input file   %
% is present.                                                            %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc;
clear;
close all;
%%
[vfile,videofolder]=uigetfile;
filename = [videofolder '/' vfile];
fullpath = videofolder;
load(filename);
%%
initial = 1; % 701 bcz we need to start from 601. 
final = length(X);
chunks = 2;%10; % 10 for grop of 15, 15-25 for grp of 30 and 60, 50 for grp of 120 or240
in= floor((final-initial)/(chunks-1));
stitch_len = 10;%30;
%%
tc=0;
for m=initial:in:final
    %disp(m);
    %lll = input('something');
    % setting initial frame for fishtracker. z is index for positions
    if m ==1
%         initial_frame=initial_frame;
        z=m;
    else
        initial_frame = (m-stitch_len)*steps;
        z=m-stitch_len;
        %lll = input('initial_frame');
    end
    
    %setting nframe 
    if m == 1
        nframe = ((m+in)*steps) + initial_frame -1 ;
    else if length(X) > (m+in)    
            nframe = (m+in)*steps ;
         else
            nframe = (length(X)-1)*steps;
         end
    end
    
    ifrm=num2str(initial_frame);    %initial frame number for naming
    nfrm=num2str(floor(nframe));    %end frame number for naming
    fish_tracker;                   %Tracking code
    data_sheet = datasheet(Q_loc_estimateX,Q_loc_estimateY,steps,frameRate);%make data sheet
    fullname = [fullpath 'tracked_' vfile(10:end-4) '_' ifrm 'to' nfrm '.mat'];
    save(fullname,'Q_loc_estimateX','Q_loc_estimateY','initial_frame','nframe','steps','data_sheet','-v7.3');

end
