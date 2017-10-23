%% find the indices with orientation more or less equal to the mode orientation
function [sameOrientIndTmp, meanYaw] = findSameOrient(orientConsider,orientThres)
% % % INPUT:
% orientConsider: yaw values in degrees at the time of transmission of each of the WIFi packets
% orientThres: threshold in degrees within which we consider orientation is more or less the same, say within 10 degrees
% % % OUTPUT:
% sameOrientIndTmp: indices of packets with ~same orientation
% meanYaw: mean of the yaw of packets with ~same orientation


%     smallOriThres = 0.005;
%     orientConsiderInd = find( abs(diff(yaw))>smallOriThres, 1, 'first' ) : find( abs(diff(yaw))>smallOriThres, 1, 'last' );
%     orientConsider = yaw(orientConsiderInd);
%     orientThres = 6; %30 ;
    histEdges = linspace(-180, 180, 1 + 360/(2*orientThres));
    [binCount,histEdges] = histcounts(orientConsider, histEdges);
    [~, maxBin] = max(binCount);
    sameOrientIndTmp = 1:length(orientConsider);find(orientConsider>histEdges(maxBin) & orientConsider<histEdges(maxBin+1));
    meanYaw = mean(orientConsider(sameOrientIndTmp));
    % sameOrientInd = dispExpDataIdTot((sameOrientIndTmp));
%     sprintf('Ratio of points with more or less same orientation is %f', length(sameOrientIndTmp)/length(orientConsider))
    % % making expdataId to consider points whose orientation is same as the mode
    % ExpDataIdTot = sameOrientInd(round( linspace(1,length(sameOrientInd),min(600,length(sameOrientInd))) )); %sameOrientInd(min(1, end-500):end);
    % % making hte expdataId and dataExpdataId the same
    % ExpDataIdTot = dispExpDataIdTot;
