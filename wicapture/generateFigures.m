%% genrating the WiCapture plots in the paper
clc; clear; close all;
lineWidthPlot = 1;
%% Figure 8(a), indoor office deployment experiments
loadData = 1;
if loadData % use previously calculated data
    load('cvprLos','trajErrSync_wicapture','trajErrSync_spotfi');
else % calculate the error from stored files
    expOffsetTot = 1:98;
    trajErrSync_wicapture = [];
    trajErrSync_spotfi = [];
    for iExpOffset = 1:length(expOffsetTot)
        expFileOffset = expOffsetTot(iExpOffset);
        rawData = load(sprintf('../wicaptureData/resultData/los_trajOvrWiFi%d.mat',expFileOffset));
        %rawDataSpotfi = load(sprintf('../../PinPoint123456/ovrWiFiFigs/spotfiTraj_los_%d.mat',expFileOffset));

        trajErrSync_wicapture = [trajErrSync_wicapture; rawData.trajErrSync(2:end)];
        %trajErrSync_spotfi = [trajErrSync_spotfi; rawDataSpotfi.trajErrSync(2:end)];
    end

    % making it to centimeter scale
    trajErrSync_wicapture = trajErrSync_wicapture*1e-1;
    %trajErrSync_spotfi = trajErrSync_spotfi*1e-1;
    save('cvprLos','trajErrSync_wicapture','trajErrSync_spotfi');
end

median(trajErrSync_wicapture)
%median(trajErrSync_spotfi)

figure;
h = cdfplot(trajErrSync_wicapture); hold all
set(h, 'LineWidth', lineWidthPlot, 'LineStyle', '-')
% h = cdfplot(trajErrSync_spotfi); hold all
% set(h, 'LineWidth', lineWidthPlot, 'LineStyle', '--')
hold off;
set(gca,'Xscale','log')
xlabel('Trajectory error (cm)')
ylabel('CDF')
title('')
legend('WiCapture', 'Location', 'NorthWest')
xlim([2e-2 800])
set(gca, 'XTickLabel', [0.1 1 10 100])
 

%% Outside experiments
% loadData = 0;
% if loadData
%     load('cvprOut','trajErrSync_wicapture','trajErrSync_spotfi');
% else
%     expOffsetTot = 1:89;%[37:53 55:73 75:111 113:125]; 171:268; 269:357;
%     trajErrSync_wicapture = [];
%     trajErrSync_spotfi = [];
%     for iExpOffset = 1:length(expOffsetTot)
%         expFileOffset = expOffsetTot(iExpOffset);
%         rawData = load(sprintf('../wicaptureData/resultData/outside_trajOvrWiFi%d.mat',expFileOffset));
%         %rawDataSpotfi = load(sprintf('../../PinPoint123456/ovrWiFiFigs/spotfiTraj_outside_%d.mat',expFileOffset));
% 
%         trajErrSync_wicapture = [trajErrSync_wicapture; rawData.trajErrSync(2:end)];
%         %trajErrSync_spotfi = [trajErrSync_spotfi; rawDataSpotfi.trajErrSync(2:end)];
%     end
% 
%     % making it to centimeter scale
%     trajErrSync_wicapture = trajErrSync_wicapture*1e-1;
%     %trajErrSync_spotfi = trajErrSync_spotfi*1e-1;
%     save('cvprOut','trajErrSync_wicapture','trajErrSync_spotfi');
% end
% 
% median(trajErrSync_wicapture)
% %median(trajErrSync_spotfi)
% 
% figure;
% h = cdfplot(trajErrSync_wicapture); hold all
% set(h, 'LineWidth', lineWidthPlot, 'LineStyle', '-')
% % h = cdfplot(trajErrSync_spotfi); hold all
% % set(h, 'LineWidth', lineWidthPlot, 'LineStyle', '--')
% hold off;
% set(gca,'Xscale','log')
% xlabel('Trajectory error (cm)')
% ylabel('CDF')
% title('')
% legend('WiCapture', 'SpotFi', 'Location', 'NorthWest')
% xlim([2e-2 800])
% set(gca, 'XTickLabel', [0.1 1 10 100])

