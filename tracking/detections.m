%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Amith Kumar U R                                                        %
% Description:                                                                   %
%                                                                                %
% This code takes a video file as input (.mov,.avi,.mp4), detects the            %
% position of individuals in each frame. The output is a (.mat) cell array,      %
% each element of which encapsulates X and Y coordinates of all individals       %
% detected in a frame. The method used for detection is called image subtraction.%
% For tracking, after detections loop_track.m is needs to be run.                %
%                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read video file from a folder
[vfile,videofolder]=uigetfile({'*.MOV';'*.avi';'*.mp4'},'VideoFileSelector');
videofile = [videofolder '/' vfile];
videoread=VideoReader(videofile);
nframe= round(videoread.Duration * videoread.FrameRate); %to get last frame number of/in the video
frameRate = videoread.FrameRate;

%% Create region of interest by clicking and dragging in the window that will pop up. Once selected, right click and copy positions
% and then press enter to continue

videoread.CurrentTime = videoread.Duration/2;
imgframe = readFrame(videoread);
figure, 
imshow(imgframe); 
H = imellipse;
 
pause

BW=createMask(H);
close ;


%% Initialize parameters for convolution
hsizeh = 10;  % tweek for better results.
sigma = 5; % tweek for better results. 
h = fspecial('log', hsizeh, sigma);


%% Initialize ..
steps = 1;           % determines resolution of the data (set it to 1 if fish is moving very fast, values(>3) may be difficult to track).
initial_frame = 1;	 % frame from which detection needs to be started
endframe = nframe-10; % frame till which detection needs to be done
imgsub = 3;      % this number (in seconds) will be used to subtract frames.
X = cell(floor(endframe-initial_frame+1),1); %to store X coordinates 
Y = cell(floor(endframe-initial_frame+1),1); %to store Y coordinates

%% Lets detect fish using image subtraction!! ;)
% i = 50;
c = 0; %counter
% Uncomment below if you want visualization of detections. Also uncomment at two further places below.
% nFrames = 20;
% vidObj = VideoWriter('detectedBeetles.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 10;
% open(vidObj);
for i=initial_frame:steps:endframe %% i value changed so that video starts from 3mins
    c=c+1;
    videoread.currentTime = i/videoread.FrameRate;
    imgframe = readFrame(videoread);   % Current frame number
    
    
    %% Lets pick frame for image subtraction
    % condition to choose frame for subtraction from the current frame
    if (nframe - 10) > ((videoread.CurrentTime+imgsub)*videoread.FrameRate)

        videoread.CurrentTime = videoread.CurrentTime + imgsub;
        imgframe1 = readFrame(videoread);% different frame to subtract
    else
        videoread.CurrentTime = videoread.CurrentTime - imgsub;
        imgframe1 = readFrame(videoread);% different frame to subtract
    end
    %% Image segmentation
    imgf = imgframe - imgframe1;      % Do Image subtraction
    sigleimg=im2single(imgf);         % make single
    sigleimg1=sigleimg(:,:,1);        % take out any one matrix out of 3(rgb)
    sigleimg1(sigleimg1 <  0.1) = 200;% filtering
    sigleimg1(BW < 1) = 200;        % mask only ROI
    
    blob_img = conv2(sigleimg1,h,'same');%Convolution to get blobs
    
    blob_img(blob_img < 0.05) = nan ;% filter blob image
%     figure, imshow(blob_img)
%     pause;
    % Check for at least one maxima other wise we assume all blobs are stationary
    if max(max(blob_img)) > 0.0001
        
        [zmax,imax,zmin,imin] = extrema2(blob_img); % find extrema
        [X{c},Y{c}] = ind2sub(size(blob_img),imax); % Store max values
        count(c,1)=length(X{c});    % number of blobs detected
        
    else
        if c == 1
            X{c} = nan;
            Y{c} = nan;
            count(c,1)=0;
        else %if no blobs detected we assume all blobs are stationary
            X{c} = X{c-1};
            Y{c} = Y{c-1};
            count(c,1)=length(X{c});
        end
    end
    count(c,2) = videoread.CurrentTime;
    cnt=count(c,1)      % Display Count
    % ------------------------------------------------------------------------
    %% Visualization (Uncomment block of code below once)
    labels=0;
%     imshow(imgframe1);
%     hold on;
%     plot(Y{c},X{c},'.','MarkerSize',10)
%     for j = 1:length(X{c})
% %         labels = labels + 1;
% %         labels=cellstr( num2str(j));
%         plot(Y{c}(j),X{c}(j),'r.','MarkerSize',10)
% %         text(Y{c}(j),X{c}(j),labels,'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize',.1)
%         
%     end
%     f=getframe(gca);

%     pause(0.01);
    

%   writeVideo(vidObj, f);
    
completed = 100*i/(endframe) %(displays in percentage)
%     completed = videoread.CurrentTime %(displays in time)

end

% close(vidObj); %uncomment for visualization

%% saving into the location where your video is present
fulname = [videoread.Path '/detected_' videoread.Name(1:end-4) '.mat'];
path = videoread.Path;
Name = videoread.Name(1:end-4) ;
save(fulname, 'X','Y','count','steps','initial_frame','endframe','path','Name','BW','frameRate');


