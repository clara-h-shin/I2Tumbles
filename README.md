## I2Tumbles (SimpleTracking V 2.0)
I2Tumbles converts bacterial movement images into quantitative tracking data, analyzes bacterial motility, 
and provides metrics such as tumbles per second to help researchers compare experimental outcomes.

This script was initially developed to support research in the Baylink Lab, but it can also be applied to a wide range of chemotaxis studies.

To use this program: 
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


#### Baylink Lab
www.baylink-lab.com
