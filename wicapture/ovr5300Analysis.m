clear; close all; clc;

expOffsetTot = [171:193 231:268];%[37:53 55:73 75:111 113:125]; 171:268; 269:357;
trajErrSync = [];
dispErrSync = [];
trajErrSyncEachData = [];
rmse = [];
rmseNormalized = [];
dispErrSyncNormalized = [];
trajErrSync_kalman = [];
dispErrSync_kalman = [];
dispErrSyncNormalized_kalman = [];
rmse_kalman = [];
for iExpOffset = 1:length(expOffsetTot)
    expFileOffset = expOffsetTot(iExpOffset);
    %rawData = load(sprintf('../../PinPoint123456/ovrWiFiFigs/pkd75_run1_trajOvrWiFi%d.mat',expFileOffset));
    rawData = load(sprintf('../../PinPoint123456/ovrWiFiFigs/ignoreCfo_trajOvrWiFi%d.mat',expFileOffset));

%     rawData = load(sprintf('../../PinPoint123456/ovrWiFiFigs/elevAp_allCsi_TrajOvrWiFi%d.mat',expFileOffset));
%     rawData = load(sprintf('../../PinPoint123456/ovrWiFiFigs/elevApTrajOvrWiFi%d.mat',expFileOffset));
%     saveas(gcf,sprintf('../../PinPoint123456/ovrWiFiFigs/pic_elevApTrajOvrWiFi%d.jpg',expFileOffset))
%     rawData = load(sprintf('../../PinPoint123456/ovrWiFiFigs/elevApTraj_aoa_OvrWiFi%d.mat',expFileOffset));
%     saveas(gcf,sprintf('../../PinPoint123456/ovrWiFiFigs/pic_elevApTraj_allCsi_OvrWiFi%d.jpg',expFileOffset))
    trajErrSync = [trajErrSync; rawData.trajErrSync];
    dispErrSync = [dispErrSync; rawData.dispErrSync];
    dispErrSyncNormalized = [dispErrSyncNormalized; rawData.dispErrSyncNormalized];

    rmse = [rmse; rawData.rmse];
    rmseNormalized = [rmseNormalized; rawData.rmseNormalized];
    
%     trajErrSyncEachData = [trajErrSyncEachData rawData.trajErrSync];
    
    sprintf('For exp %d, normalized rmse is %f', expFileOffset, rawData.rmseNormalized)
    


%     accInd = 11:13;
%     timeInd = 1;
%     accData = rawData.ovrTmp(rawData.nearTimeInd, accInd);
%     accData(:,2) = accData(:,2) - 9.9; % removing g value
%     timestepData_s = diff(rawData.ovrTmp(rawData.nearTimeInd, timeInd))*1e-6;
%     disp_m = diff(rawData.trajEstShift)*1e-3;
%     trueTraj_m = rawData.trajOvrSync*1e-3;
% %     [traj_kalmanDisp_m, MM] = kalmanDispDiffAcc_3D(disp_m, accData(1:end-1,:), timestepData_s, trueTraj_m);
% %     traj_kalmanDisp = traj_kalmanDisp_m*1e3;
%     [traj_kalmanDisp, MM] = kalmanDispDiffAlone_3D(diff(rawData.trajEstShift), rawData.trajOvrSync);
%     f = figure;
%     plot3( traj_kalmanDisp(:,1), traj_kalmanDisp(:,2), traj_kalmanDisp(:,3), 'd-', 'LineWidth', 2); hold on; 
%     plot3( rawData.trajOvrSync(:,1), rawData.trajOvrSync(:,2), rawData.trajOvrSync(:,3), '-', 'LineWidth', 4);  hold off; axis equal; axis tight; 
%     legend('from WiFi', 'true path', 'Location','NorthOutside','Orientation','horizontal')
%     xlabel('x(mm)'); ylabel('y(mm)'); zlabel('z(mm)')
%     h = get(gca,'DataAspectRatio');
%     if h(3)==1
%           set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
%     else
%           set(gca,'DataAspectRatio',[1 1 h(3)])
%     end
%     trajErrSync_kalman = [trajErrSync_kalman; norms(traj_kalmanDisp - rawData.trajOvrSync, 2, 2)];
%     dispErrSync_kalman = [dispErrSync_kalman; norms(diff(traj_kalmanDisp) - diff(rawData.trajOvrSync), 2, 2)];
%     dispErrSyncNormalized_kalman = [dispErrSyncNormalized_kalman; norms(diff(traj_kalmanDisp) - diff(rawData.trajOvrSync), 2, 2)./norms(diff(rawData.trajOvrSync), 2, 2)];
%     rmse_kalman = [rmse_kalman; sqrt(mean(trajErrSync_kalman.^2))];
end
median(trajErrSync)
median(rmse)
median(rmseNormalized)
save(sprintf('../../PinPoint123456/cvprResults/trajErrIgnoreCfo.mat'))

%%
pktConsider = 1:5:size(rawData.ovrTmp,1);
accDataRaw = rawData.ovrTmp(pktConsider, accInd);
accData = accDataRaw;
accData(:,2) = accData(:,2) - 9.9; % removing g value
timestepData_s = diff(rawData.ovrTmp(pktConsider, timeInd))*1e-6;
timestepMat = repmat(timestepData_s,1,3);
v = cumsum(accData(1:end-1,:).*timestepMat);
s = v.*timestepMat + 0.5*accData(1:end-1,:).*(timestepMat.^2);
figure; plot(s)
title('s')

figure; plot(trueTraj_m); title('true s')
%%
sprintf('75th percentile of normalized displacement error %f', prctile(dispErrSyncNormalized,75))
sprintf('75th percentile of normalized trajectory error %f', prctile(trajErrSync,75))
sprintf('75th percentile of displacement estimation error %f', prctile(dispErrSync,75))
sprintf('rmse normalized is %f', rmseNormalized)
sprintf('rmse is %f', rmse)

figure; 
cdfplot(rmse); xlabel('RMSE (mm)')
figure;
cdfplot(rmseNormalized); xlabel('RMSE (mm) normalized to path length')
figure;
cdfplot(dispErrSync); xlabel('displacement estimation error (mm)')
set(gca,'XScale','log')
figure;
cdfplot(dispErrSyncNormalized); xlabel('displacement estimation error normalized to actual displacement (mm)')
set(gca,'XScale','log')

figure; 
cdfplot(rmse_kalman); xlabel('kalman RMSE (mm)')
figure;
cdfplot(dispErrSync_kalman); xlabel('kalman displacement estimation error (mm)')
set(gca,'XScale','log')
figure;
cdfplot(dispErrSyncNormalized_kalman); xlabel('kalman displacement estimation error normalized to actual displacement (mm)')
set(gca,'XScale','log')

sprintf('50th percentile of normalized displacement error %f', prctile(dispErrSyncNormalized,50))
sprintf('50th percentile of normalized trajectory error %f', prctile(trajErrSync,50))
sprintf('50th percentile of displacement estimation error %f', prctile(dispErrSync,50))
sprintf('median RMSE is %f', median(rmse))

sprintf('kalman 50th percentile of normalized displacement error %f', prctile(dispErrSyncNormalized_kalman,50))
sprintf('kalman 50th percentile of normalized trajectory error %f', prctile(trajErrSync_kalman,50))
sprintf('kalman 50th percentile of displacement estimation error %f', prctile(dispErrSync_kalman,50))
sprintf('kalman median RMSE is %f', median(rmse_kalman))

%% code for finding displacement for non-consecutive packets
addpath('./automowfilter');
succPktDiffConsiderTot = 1;

rot2d = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];

dispErrSyncTot = [];
dispMagTot = [];
for iExpOffset = 1:length(expOffsetTot)
    expFileOffset = expOffsetTot(iExpOffset);
    aTmp = load(sprintf('../../PinPoint123456/ovrWiFiFigs/trajOvrWiFi%d.mat',expFileOffset));
%     close;

    for succPktDiff = succPktDiffConsiderTot
        for iDiffStart = 1:succPktDiff
            dispErrSuccDiff = norms(diff(aTmp.trajEstShift(iDiffStart:succPktDiff:end,:)) - diff(aTmp.trajOvrSync(iDiffStart:succPktDiff:end,:)), 2, 2);
            dispSuccDiff = norms(diff(aTmp.trajOvrSync(iDiffStart:succPktDiff:end,:)), 2, 2);

            dispErrSyncTot = [dispErrSyncTot; dispErrSuccDiff];

            % dispMag = norms(diff(aTmp.trajOvrSync), 2, 2);
            dispMagTot = [dispMagTot; dispSuccDiff];
        end
    end
end
median(dispErrSyncTot)
figure(3); cdfplot(dispErrSyncTot)

%%
edges = (1:2:50);
dispMagDiscrete = discretizeMani(dispMagTot,edges);
dispMagDiscrete = edges(dispMagDiscrete+1);
figure(2); hist(dispMagTot, edges)
figure(4); boxplot(dispErrSyncTot, dispMagDiscrete, 'Whisker', Inf)


%% old code for orientation correction
% tilefigs
close all;
for iExpOffset = 1:length(expOffsetTot)
    expFileOffset = expOffsetTot(iExpOffset);
    aTmp = load(sprintf('../../PinPoint123456/ovrWiFiFigs/trajOvrWiFi%d.mat',expFileOffset));
    
    aTmp.trajEstSync;
    aTmp.trajOvrSync;
    
    oriInd = 7:10;
    [pitch_rad, yaw_rad, roll_rad] = quat2angle(aTmp.ovrTmp(:,oriInd));
    yaw = rad2deg(yaw_rad); 
    pitch = rad2deg(pitch_rad);
    roll = rad2deg(roll_rad);
    figure(2); plot(aTmp.trajEstSync(:,1), aTmp.trajEstSync(:,2), aTmp.trajOvrSync(:,1), aTmp.trajOvrSync(:,2) ); axis equal;
    
%     trajEst_yawCorrect = rot2d(-mean(yaw_rad))*aTmp.trajEstSync.';
%     trajEst_yawCorrect = trajEst_yawCorrect.';
    trajEst_yawCorrect = aTmp.trajEstSync;
    for iPkt = 1:size(aTmp.trajEstSync,1);
        trajEst_yawCorrect(iPkt,:) = aTmp.trajEstSync(iPkt,:)*rot2d(yaw_rad(iPkt));
    end
    figure(3); plot(yaw);
    
    figure(4); plot(trajEst_yawCorrect(:,1), trajEst_yawCorrect(:,2), aTmp.trajOvrSync(:,1), aTmp.trajOvrSync(:,2) ); axis equal;
    title('yaw corrected')
    
    f2 = figure(5);
    plot( aTmp.trajEstShift(:,2), -aTmp.trajEstShift(:,1), 'd-', 'LineWidth', 2); hold on; 
    plot( aTmp.trajOvrSync(:,2), -aTmp.trajOvrSync(:,1), '-', 'LineWidth', 4);  hold off; axis equal; axis tight; 

end
