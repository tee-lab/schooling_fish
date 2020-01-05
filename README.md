![Etroplus suratensis](https://teelabiisc.files.wordpress.com/2019/03/fish-4.jpg)
# schooling_fish

## This repository contains codes and data used for analysis in the article "Noise-Induced Schooling of Fish". The first release is archived at [![DOI](https://zenodo.org/badge/159283981.svg)](https://zenodo.org/badge/latestdoi/159283981)



**Note: Codes have been tested in a Linux and Windows 10 machine. A few adjustments may be required for Mac users. For all the inputs set the working directory in MATLAB to where the input files are saved.**

# Codes for detection and tracking

1. **Detection:** Detection of fish from a video can be done using the *MATLAB* code __detection.m__ in ``/tracking/detections.m``.
    1. Run the file __detections.m__ which will prompt you to select the video for detection.
    2. Use sample video __sample_vid.avi__ in the folder ``/sample`` as input for the code.
    3. A figure pops up displaying a frame of the video. Click and drag to select the region where you want the fish to be detected, right click and select ``copy position``, minimize the window and hit ENTER on the *MATLAB* console to begin tracking. 
    >*-- Note: It could be the entire arena where the fish swim or a small region where they are isolated or is of interest, also make sure not to close the window after selecting the region, just minimize it.
    
    -- The time to run this code would depend on the length of the video and your system. For the sample video provided by us it would take approximately 40 mins.*
    
    4. Individual fish are detected and their coordinates at every frame is stored in the folder ``/sample`` as *detected_sample_vid.mat* and *detected_sample_vid.csv*.
     >Note: Change 'steps' in the code to steps > 1 for skipping frames between detection.
 
2. **Tracking:** A *MATLAB* code ``/tracking/loop_track.m`` can be used for tracking individuals between frames. This code uses the Kalman-filter code which we downloaded from [Mr Student Dave’s website tutorial](http://studentdavestutorials.weebly.com/). Compared to ``/tracking/detections.m`` this code can be expected to take substantially less time.
    1. Use *detected_sample_vid.mat* as the input when prompted.
    2. The output is a set of files with name as in *tracked_sample_vid_###to###.mat*, where *###* is the start-frame number and end-frame number respectively in the .mat file. These files are also stored in the folder ``/sample``.

3. **Stitching different output files:** A *MATLAB* code  ``/tracking/meta_analysis.m`` can be used to stitch the output data from above code __loop_track.m__ into one single file called *rawdata_track.csv*. This output will be stored in the folder ``/sample``.
   >*Note: Owing to size of video files, we have provided only a sample video for group size (N) 60.*

# Codes for Figure 2

1. **Generate polarization time-series for sample video:** The *MATLAB* code ``/data_analysis/op_calculate.m`` can be used to compute the polarization/order parameter.
    1. The input for this code is *rawdata_track.csv* in the ``/sample``. For the data used in the paper select the ones provided in ``/data/raw_data/``. 
    For e.g,``/data/raw_data/15/raw_data15tr1.csv`` as input file.
    3. The output vel_x, vel_y are saved as *<file_name>_vel_x_y.mat* for further analysis by the users.
    
**To generate Fig 2**, we used different N and multiple trials for each N. To reproduce these plots, one must first generate the corresponding complete-time series (i.e. each trial of around 50 minutes at ~8 fps). To do this, run the rawdata files ``/data/raw_data`` for different N. 
The time series and histogram plots are generated using vel_x and vel_y for various group size and trials.

# Codes for Figure 3

1. Drift and Diffusion functions can be computed by running the *MATLAB* code ``/data_analysis/vector_drift_diff.m``

Input required: (vel_x.mat and vel_y.mat) computed from vel_x, vel_y computed as explained in section ‘Codes for Figure 2’ or synthetically generated using ``/simulations/interact.m`` as explained in next section, ‘Codes for Figure 4’.

# Codes for Figure 4

1. A *MATLAB* code ``/simulations/interact.m`` used for simulating the pairwise and ternary interaction models.
    1. For pairwise set the parameter controlling ternary interaction rate (h) = 0 in the code with s = 0.25 and p = 4 in line 11 and sigma=3 in line 12.
        1. The outputs (polarization vectors) vel_x.mat and vel_y.mat are saved when required.
        2. The saved output can be used to plot drift and diffusion functions using the codes in ``/data_analysis/vector_drift_diff.m``
    2. For Ternary set the parameter controlling ternary rate (h)=0.3 in the code with s = 0.25 and p = 0.01 in line 11 and sigma=3 in line 12.(line 11).
        1. The outputs (polarization vectors) vel_x.mat and vel_y.mat are saved when required. 
        2. The saved output can be used to plot drift and diffusion functions using the codes in ``/data_analysis/vector_drift_diff.m``

# GA simulations

1. In the folder ``/Simulations/GA`` you will find the codes to run Genetic algorithm optimization for identifying rates of interaction between fish and the order of alignment interactions.
    1. In the sub-folder, ``/Copyinginteractions`` the alignment model is based on a hierarchy of models where the fish copy the **direction from another fish** from a set k other fish.
    2. In the sub-folder, ``/averaginginteractions`` the alignment model is based on a hierarchy of models where the fish choose the **average direction of k other fish**.
    
    
For both the cases, Run the __MainFile.m__ to begin the optimization.
For standalone codes of just the model (copying and averaging type), refer to the folder ``/GA/standalone``.
 
 
 
#### Note: codes have been tested in Linux and Windows 10 machine. A few adjustments may be required for Linux and Mac.
