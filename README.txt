Download and extract web.stanford.edu/~mkotaru/wicaptureData.zip file into this directory.

./wicapture/generateFigures.m plots the CDF plot of error for trajectories estimated by WiCapture system for experiments conducted in indoor office deployment illustrated in Figure 8(a) in the paper http://openaccess.thecvf.com/content_cvpr_2017/papers/Kotaru_Position_Tracking_for_CVPR_2017_paper.pdf. 

./wicapture/genTrajectoryLos.m generates the trajectories estimated by WiCapture, the data that is plotted in the above CDF plot. The code is tested in MATLAB R2014b.

./supplementaryMaterial.pdf provides a detailed description of how to use multiple receive antennas and multiple frequencies to obtain higher accuracy in motion tracking and also describes a one-time calibration procedure to compensate for transceivers’ chain responses.