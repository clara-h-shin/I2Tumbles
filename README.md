## I2Tumbles (SimpleTracking V 2.0)
I2Tumbles converts bacterial movement images into quantitative tracking data, analyzes bacterial motility, 
and provides metrics such as tumbles per second to help researchers compare experimental outcomes.

This script was initially developed to support research in the Baylink Lab, but it can also be applied to a wide range of chemotaxis studies.

###To use this program:### 
1. Run *SimpleTracking.m*
  - This script converts bacterial movement images into x and y coordinates, saved as SimpleTrackingOutput.mat.
  - You will need to have tiff images converted from a video to run this script properly. To learn more about converting the images, please refer to *swim_tracking_tutorial_by_kailie.docx*.
  - Set the parameters before running. For more detailed instructions, refer to *SimpleTracking_Instructions.pdf*.
2. Run *CellTrackAnalysis.m*
  - This script reads in SimpleTrackingOutput.mat or any other mat file from TrackingGUI, and provides summarized and individual metrics of the experimental group.
  - Set the parameters before running. For more detailed instructions, refer to *CellTrackAnalysis_Instructions.pdf*.
3. [Optional] Run *PlotTrajectory.m*
  - This script can be used for checking the validity of the outcomes from *PlotTrajectory.m*.
  - Set the parameters before running. For more detailed instructions, refer to *PlotTrajectory_Instructions.pdf*.

###Inputs and Outputs###
1. *SimpleTracking.m*
   - Input: .tif image files converted by Fiji/ ImageJ from a video. Should be more than 1 image.
   - Outputs:
     - Prints out the number of tracks
     - SimpleTrackingoutput.mat: A mat file with metadata to be used in CellTrackAnalysis.m for tumble analysis. This can be replaced by the output .mat file from TrackingGUI
     - tracked_positions.xlsx: An excel file with X, Y coordinates of bateria, metrics including instantaneous speed.
     - Visualization: The first image with bacterium movement tracking drawn on
2. *CellTrackAnalysis.m*
   - Input: SimpleTrackingoutput.mat or .mat file from TrackingGUI
   - Outputs:
     - Prints out the summarized (mean and median) metrics (i.e. Number of Tumbles, Tumbles Per Second, Tumble Angle, etc.)
     - CellTrackAnalysis_results.xlsx: An excel file with metrics per bacterium and summarized metrics
     - Visualization: Histograms of the metrics and a scatterplot of Tumbles Per Second vs Track Duration
3. *PlotTrajectory.m*
   - Input: SimpleTrackingoutput.mat or .mat file from TrackingGUI
   - Outputs:
     - If bacteria_idx = 0, a visualization of randomly selected 20 bacteria tracks with tumble location in blue dots
     - If bacteria_idx = bacteria id, a larger track with tumble location of a bacterium with graphs of linear and angular velocities marked with tumble location

#### Baylink Lab
www.baylink-lab.com
