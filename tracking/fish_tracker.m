%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This program will assign the tracks from the detected positions.
% Kalman filter tracking algo has been used in this code. I have modified
% the original code a bit so that works better for our setup.

%% This code is originally written by StudentDave (Donno exact name :p).
% You can get original code from below link:
% http://studentdavestutorials.weebly.com/kalman-filter-with-matlab-code.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define main variables for KALMAN FILTER! :P
dt = 1;  %our sampling rate
z = 1;
S_frame = z;%floor(d); % find(cellfun(@length, X)>11,1); %starting frame

%now, since we have multiple flies, we need a way to deal with a changing
%number of estimates! this way seems more clear for a tutorial I think, but
%there is probably a much more efficient way to do it.

u = 0; % define acceleration magnitude to start
HexAccel_noise_mag = 3; %process noise: the variability in how fast the Hexbug is speeding up (stdv of acceleration: meters/sec^2)
tkn_x = 1;  %measurement noise in the horizontal direction (x axis).
tkn_y = 1;  %measurement noise in the horizontal direction (y axis).
Ez = [tkn_x 0; 0 tkn_y];
Ex = [dt^4/4 0 dt^3/2 0; ...
    0 dt^4/4 0 dt^3/2; ...
    dt^3/2 0 dt^2 0; ...
    0 dt^3/2 0 dt^2].*HexAccel_noise_mag^2; % Ex convert the process noise (stdv) into covariance matrix
P = Ex; % estimate of initial Hexbug position variance (covariance matrix)

%% Define update equations in 2-D! (Coefficent matrices): A physics based model for where we expect the HEXBUG to be [state transition (state + velocity)] + [input control (acceleration)]
A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]; %state update matrice
B = [(dt^2/2); (dt^2/2); dt; dt];
C = [1 0 0 0; 0 1 0 0];  %this is our measurement function C, that we apply to the state estimate Q to get our expect next/new measurement




%% initize result variables
Q_loc_meas = []; % the fly detecions  extracted by the detection algo
%% initize estimation variables for two dimensions
Q= [X{S_frame} Y{S_frame} zeros(length(X{S_frame}),1) zeros(length(X{S_frame}),1)]';
Q_estimate = nan(4,25000);
nbrof_frames= floor((nframe-initial_frame)/steps) +1;
Q_estimate(:,1:size(Q,2)) = Q;  %estimate of initial location estimation of where the flies are(what we are updating)
Q_loc_estimateY = nan(nbrof_frames,5000); %  position estimate
Q_loc_estimateX= nan(nbrof_frames,5000); %  position estimate
P_estimate = P;  %covariance estimator
strk_trks = zeros(1,100000);  %counter of how many strikes a track has gotten
nD = size(X{S_frame},1); %initize number of detections
nF =  find(isnan(Q_estimate(1,:))==1,1)-1 ; %initize number of track estimates
%for each frame

c=0;%floor(d-1);
%for t = S_frame:((length(f_list)))-1 
%  z=initial_frame;%c=c+1;
%while ~isDone(videoread);
%%
% nFrames = 20;
% vidObj = VideoWriter('videoTracked.avi');
% vidObj.Quality = 100;
% vidObj.FrameRate = 20;
% open(vidObj);
nframe = endframe;
z=z-1;
 for i=initial_frame:steps:nframe
     c=c+1;
   z=z+1;
  % c=c+1;
    
    % load the image
   % img_tmp = double(imread(f_list(t).name));
   
    % make the given detections matrix
    Q_loc_meas = [X{z} Y{z}];
    
    
    %% do the kalman filter
    % Predict next state of the flies with the last state and predicted
    % motion.
    nD = size(X{z},1); %set new number of detections
    for F = 1:nF
         Q_estimate(:,F) = A * Q_estimate(:,F) + B * u;
        
    end
    
    %predict next covariance
    P = A * P* A' + Ex;
    % Kalman Gain
    K = P*C'*inv(C*P*C'+Ez);
    
    
    %% now we assign the detections to estimated track positions
    %make the distance (cost) matrice between all pairs rows = tracks, coln
    %= detections
    est_dist = pdist([Q_estimate(1:2,1:nF)'; Q_loc_meas]);
%     est_dist = est_distn(1,~isnan(est_distn));
    est_dist = squareform(est_dist); %make square
    est_dist = est_dist(1:nF,nF+1:end) ; %limit to just the tracks to detection distances 
    est=short_estdist(est_dist,nF);
    set(0,'RecursionLimit',1000);
    [asgnt, cost] = munkres(est);  % assignmentoptimal(est); %do the assignment with hungarian algo
    asgnt = asgnt'; 
    clear indx;
    indx=find(~isnan(est_dist(:,1)));
    asgn = zeros(1,nF);
    asgn(indx)=asgnt; 
%     asgn = zeros(1,nF);
%     asgn(indx)=asgnt; 
    % ok, now we check for tough situations and if it's tough, just go with
    % estimate and ignore the data
    %make asgn = 0 for that tracking element
    
    %check 1: is the detection far from the observation? if so, reject it.
    rej = [];
    for F = 1:nF
        if asgn(F) > 0
            rej(F) =  est_dist(F,asgn(F)) <  15; %set aproximate value such that path shouldn't overshoot.
        else
            rej(F) = 0;
        end
    end
    asgn = asgn.*rej;
        
    %apply the assingment to the update
    k = 1;
    for F = 1:length(asgn)
        if asgn(F) > 0
            Q_estimate(:,F) = Q_estimate(:,F) + K * (Q_loc_meas(asgn(F),:)' - C * Q_estimate(:,F));
        end
        k = k + 1;
    end
    
    % update covariance estimation.
    P =  (eye(4)-K*C)*P;
    
    %% Store data
    Q_loc_estimateX(c,1:nF) = Q_estimate(1,1:nF);
    Q_loc_estimateY(c,1:nF) = Q_estimate(2,1:nF);
    
    %ok, now that we have our assignments and updates, lets find the new
    %detections and lost trackings
    
    %find the new detections. basically, anything that doesn't get assigned
    %is a new tracking
    new_trk = [];
    new_trk = Q_loc_meas(~ismember(1:size(Q_loc_meas,1),asgn),:)';
    if ~isempty(new_trk)
        Q_estimate(:,nF+1:nF+size(new_trk,2))=  [new_trk; zeros(2,size(new_trk,2))];
        nF = nF + size(new_trk,2);  % number of track estimates with new ones included
    end
    
    %give a strike to any tracking that didn't get matched up to a
    %detection
    no_trk_list =  find(asgn==0);
    if ~isempty(no_trk_list)
        strk_trks(no_trk_list) = strk_trks(no_trk_list) + 1;
    end
    
    %if a track has a strike greater than 6, delete the tracking. i.e. make
    %it nan first vid = 3
    bad_trks = find(strk_trks > 2);
    Q_estimate(:,bad_trks) = NaN;
%--------------------------------------------------------------------------    
%%     % Just for Visualization
%     videoread.CurrentTime = i/frameRate;
%     img_tmp=readFrame(videoread);
% %     img_tmpd=double(img_tmp);
% %     img = img_tmpd(:,:,1);
% %    clf
%     %img = imread(f_list(t).name);
%     imshow(img_tmp);
%     hold on;
%    plot(Y{c}(:),X{c}(:),'r.'); % the actual tracking
%    % T = size(Q_loc_estimateX,2);
%     Ms = [1]; %marker sizes
%     c_list = ['r' 'b' 'g' 'c' 'm' 'y'];
%     for Dc = 1:nF
%         if ~isnan(Q_loc_estimateX(z,Dc))
%             Sz = mod(Dc,1)+1; %pick marker size
%             Cz = mod(Dc,6)+1; %pick color
%             if c < 21
%                 st = c-1;
%             else
%                 st = 19;
%             end
%             tmX = Q_loc_estimateX(1:z,Dc);%(c-st:c,Dc);
%             tmY = Q_loc_estimateY(1:z ,Dc);%(c-st:c,Dc);
%             plot(tmY,tmX,'.-','markersize',Ms(Sz),'color',c_list(Cz),'linewidth',1)
%             %plot(centrd(c,2),centrd(c,1),'b.')
%             axis off
%         end
%     end
%     f=getframe(gca);
%     pause(0.01);
%     writeVideo(vidObj, f);
%     
%     %}
completed = 100 * c/(((nframe-initial_frame)/steps) + 1)
 end
%  close(vidObj);
 Q_loc_estimateX(Q_loc_estimateX == 0)=nan;
 Q_loc_estimateY(Q_loc_estimateY == 0)=nan;