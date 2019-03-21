Files that need to be run in the following sequence are:
1) detections.m
2) loop_track.m
3) meta_analysis.m

The rest are dependencies. 

detection.m uses image subtraction method for detecting individuals. Specific parameters that need to be tweaked for best 
results are highlighted within the codes using inline commenting. Image subtraction method is ideal for situations where 
background is stationary with respect to individuals and individuals are constantly moving such as fish. A sufficient contrast
between the background and object is preferable. In case of heterogeneous but stationary background one may decrease the 
filtering to get blobs of detection (Line no. 77).

loop_track.m divided the input from detections.m (detected*.mat) into chunks and runs the kalman-filter algorithm code fish_tracker.m and 
generates trac*.mat files to be stiched later together by meta_analysis.m. The final output is trac_metadat.mat

# Important note: The loop_track.m is designed for analyzing large videos. For this particular sample_vid.avi it gives error while saving the last file which can be ignored and the remaining data may be stitched using meta_analysis.m. This is because
the length of the video is small.
