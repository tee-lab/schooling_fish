# schooling_fish

## This repository contains codes and data for "Finite-size noise induced schoolong in fish". It is divided into three folders, tracking, data_analysis and simulations. 

# Codes used for Figure 2

  1. A matlab code (/tracking/detections.m) used for detecting fish. Consider file (/sample/sample_vid.avi). Use it as input for the code. The output which is coordinates for all detected individuals for every frame will be stored back in the folder (/sample) as detected_sample_vid.mat.
  
  2. A matlab code (/tracking/loop_track.m) used for tracking individuals between frames. The input for this code is the output file (detected_sample_vid.mat) generated using the above code (/tracking/detections.m). The output is a collection of files (tracked_sample_vid_*to*.mat), where * indicates starting frame number to end frame number in that file. These files will also be stored in the folder (/sample).
  
  3. A matlab code (/tracking/meta_analysis.m) used for stiching the output data from above code (loop_track.m) into one single file called track_metadata.mat. This output will be stored in the folder (/sample)
  
  4. A matlab code (/data_analysis/op_calculate.m) used for calculating polarization/order parameter. The input for this code is the output (trac_metadata.m) generated using the above code (meta_analysis.m). The output (vel_x, vel_y, op) may be saved as .mat by the user for further analysis.
  
  
# Codes used for Figure 3

  1. A matlab code (/data_analysis/vector_drift_diff.m) used for calculating drift and diffusion functions of the vectors (vel_x.mat and vel_y.mat) that were calculated using the previous code (op_calculate.m). 

