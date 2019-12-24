# schooling_fish

## This repository contains codes and data used for analysis in the article "Noise-Induced Schooling of Fish". 


Note: codes have been tested in a Linux and Windows 10 machine. A few adjustments may be for using in Mac. 

# Codes for detection and tracking

1. **Detection:** Detection of fish from a video can be done using the matlab code (/tracking/detections.m). 
    1. Run the file ‘detections.m’ which will prompt you to select the video for detection.
    2. Use sample video ‘sample_vid.avi’ in the folder ‘/sample’ as input for the code.
    3. A figure pops up with the first frame of the video. Select the region where you want the fish to be detected and then press ENTER. (Note: it could be the entire arena where the fish swim or a small region where they are isolated or is of interest).
    4. Individual fish are detected and their coordinates at every frame are stored in the folder (/sample) as ‘detected_sample_vid.mat’.
2. **Tracking:** A matlab code (/tracking/loop_track.m) can be used for tracking individuals between frames. This code uses the Kalman-filter code which we downloaded from Mr Student Dave’s website tutorial (http://studentdavestutorials.weebly.com/).
    1. Use detected_sample_vid.mat as the input when prompted.
    2. The output is a set of files with name as in (tracked_sample_vid_###to###.mat), where ### are the numbers corresponding to the start-frame number and end-frame number respectively in the .mat file. These files are also stored in the folder (/sample).
3. **Stitching different output files:** A matlab code (/tracking/meta_analysis.m) can be used to stitch the output data from above code (loop_track.m) into one single file called track_rawdata.mat. This output will be stored in the folder (/sample).

Owing to size of video files, we have provided only a sample video for group size (N) 60.

# Codes for Figure 2

1. **Generate polarisation time-series for sample video:** A matlab code (/data_analysis/op_calculate.m) can be used to compute the polarization/order parameter.
    1. The input for this code is (trac_rawdata.m) in the folder/sample.
    2. The output (vel_x, vel_y, op) are saved as .mat by the user for further analysis.
    
**To generate Fig 2**, we used different N and multiple trials for each N. To reproduce these plots, one must first generate the corresponding complete-time series (i.e. each trial of around 50 minutes at ~8 fps). To do this, run the rawdata files (/data/raw_data) for different N. 
The time series and histogram plots are generated using vel_x.mat and vel_y.mat for various group size and trials.

# Codes for Figure 3

1. Drift and Diffusion functions can be computed by running the matlab code (/data_analysis/vector_drift_diff.m)

Input required: (vel_x.mat and vel_y.mat) computed from vel_x, vel_y computed as explained in section ‘Codes for Figure 2’ or synthetically generated using (/simulations/interact.m) as explained in section ‘Codes for Figure 4’.

# Codes for Figure 4

1. A matlab code (/simulations/interact.m) used for simulating the pairwise and ternary interaction models.
    1. For pairwise set the parameter controlling ternary interaction rate (h) = 0 in the code with s = 0.25 and p = 4 in line 11 and sigma=3 in line 12.
        1. The outputs (polarization vectors) vel_x.mat and vel_y.mat are saved when required.
        2. The saved output can be used to plot drift and diffusion functions using the codes in ‘/data_analysis/vector_drift_diff.m’.
    2. For Ternary set the parameter controlling ternary rate (h)=0.3 in the code with s = 0.25 and p = 0.01 in line 11 and sigma=3 in line 12.(line 11).
        1. The outputs (polarization vectors) vel_x.mat and vel_y.mat are saved when required. 
        2. The saved output can be used to plot drift and diffusion functions using the codes in ‘/data_analysis/vector_drift_diff.m’.

# GA simulations

In the folder /Simulations/GA you will find the codes to run Genetic algorithm optimization for identifying rates of interaction between fish and the order of alignment interactions.
    1. In the sub-folder, /Copyinginteractions the alignment model is based on a hierarchy of models where the fish copy the direction from another fish from a set k other fish.
    2. In the sub-folder, /averaginginteractions the alignment model is based on a hierarchy of models where the fish choose the average direction of k other fish.
    
    
For both the cases, Run the mainfile.m to begin the optimization.
For standalone codes of just the model (copying and averaging type), refer to the folder /GA/standalone.
 
 
 
Note: codes have been tested in a Windows 10 machine. A few adjustments may be required for Linux and Mac.




  

