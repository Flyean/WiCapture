%% generate WiFi trajectories for each of 98 experiments conducted in 'indoor office deployment' scenario
clear; clc; close all;

for iExpOffset = 1:89
load(sprintf('../wicaptureData/outsideData/outside%d.mat',iExpOffset), 'trajOvrSync', 'csi_trace_all', 'csi_trace_disp_all', 'nSegments', 'segStart', 'segEnd', ...
    'nSegmentsDisp','ExpDataIdTot','dispExpDataIdTot','yaw','yaw_ExpDataId')

% apply WiCapture algorithm on CSI to obtain the WiFi trajectory estimates
trajEst = wicaptureAlgo(csi_trace_all, csi_trace_disp_all, segStart, segEnd,nSegmentsDisp,ExpDataIdTot,dispExpDataIdTot,yaw,yaw_ExpDataId);
trajEst = trajEst(:,1:2); % each row of trajEst is the estimate of WiFi position for the particular packet

% WiFi transmitter is rigidly mounted on Oculus headset but their coordinate systems are not aligned. So, the WiFi position etomates are translated and rotated, but not scaled, to match the trajectory estimated by Oculus.
trajOvrOrig = trajOvrSync;
trajOvrSync(:,1) = -trajOvrOrig(:,1);
trajOvrSync(:,2) = trajOvrOrig(:,3);
trajOvrSync = trajOvrSync(:,1:2);
trajEstSync = trajEst;
B = trajEstSync'*trajOvrSync;
[BU,~,BV] = svd(B);
rotMat = BU*BV';
trajEstShift = trajEstSync*rotMat;

% 
% f2 = figure(7);
% plot( -trajEstShift(:,1), -trajEstShift(:,2), 'd-', 'LineWidth', 2); hold on; 
% plot( -trajOvrSync(:,1), -trajOvrSync(:,2), '-', 'LineWidth', 4);  hold off; axis equal; axis tight; 
% legend('from WiFi', 'true path', 'Location','NorthOutside','Orientation','horizontal')
% xlabel('motion alog x-axis in mm'); ylabel('motion along y-axis in mm')
% % set(gca,'position',[0 0 1 1],'units','normalized')
% % saveas(f2,sprintf('../../PinPoint123456/ovrWiFiFigs/trajOvrWifi%d.png',expFileOffset))
% % saveas(f2,sprintf('/Users/mani/Dropbox/ovrTraces/trajOvrWifi%d.png',expFileOffset))

% % estimating the errors
trajErrSync = norms(trajEstShift - trajOvrSync, 2, 2);
dispErrSync = norms(diff(trajEstShift) - diff(trajOvrSync), 2, 2);
dispErrSyncNormalized = norms(diff(trajEstShift) - diff(trajOvrSync), 2, 2)./norms(diff(trajOvrSync), 2, 2);
% figure; cdfplot(dispErrSync)
prctile(dispErrSyncNormalized,75)
prctile(trajErrSync,75)
prctile(dispErrSync,75)
rmse = sqrt(mean(trajErrSync.^2))
rmseNormalized = rmse/sum(norms(diff(trajOvrSync),2,2))  %norms(trajEstShift(movInd(1):movInd(2),:) - trajOvrSync(movInd(1):movInd(2),:),2,2)/norms(diff(trajOvrSync(movInd,:)),2,2);
sprintf('total length is %f', sum(norms(diff(trajOvrSync),2,2)) )

% % saving the WiCapture positioning errors to a file
save(sprintf('../wicaptureData/resultData/outside_trajOvrWiFi%d.mat',iExpOffset))

end