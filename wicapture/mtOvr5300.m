clearvars -except wiredCalib; 
clear;
close all; clc;
warning('off','all');
% % % Experiment PHYSICAL PARAMETERS
use22 = 0;
if use22
    fc = 2.2e9;
else
    fc = 5.63e9; % 2.2e9; % 
end
M = 3;
fs = 40e6;
c = 3e8;
d = 2.6e-2;%2.4e-2;
dTx = 2.6e-2;
SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % subcarriers at which information is available
SubCarrIndPhsCor = SubCarrInd - SubCarrInd(1);
N = length(SubCarrInd);
subCarrSize = 128;
fgap = 312.5e3;
INTERP = 1;
lambda = c/fc;  
T_3Tx = 3;
tofInd = 1;
aoaInd = 2;
dodInd = 3;
alphaInd = 4;
tofTotInd = 6;

% testing if I am getting enough resolution by doing a parameter sweep
expOffsetTot = 171; %[78 81];
smoothLTot = [0.75];
nCompsTot = 3; [3 4 5 6 7];
bwTot = [ 58 ]; %[64 128 256 512];
smoothKTot = [3]; %[2 3 5 7];
nLocConsiderTot = [Inf];
regScaleTot = 1; %[0 logspace(0,1.5,5)];
layoutTot = 1; 1:40;
snrTot = Inf; [15 20 30 40];

% sweepParamArray = {smoothLTot, nCompsTot, bwTot, smoothKTot, nLocConsiderTot, regScaleTot};
% sweepParamTotCell = cell(1,numel(sweepParamArray));
% [sweepParamTotCell{:}] = ndgrid(sweepParamArray{:});
% sweepParamTot = Inf*ones(numel(sweepParamTotCell{1}), numel(sweepParamArray));
% for iSweepArray = 1:numel(sweepParamArray)
%     sweepParamTot(:, iSweepArray) = sweepParamTotCell{iSweepArray}(:); 
% end

bestCompDec = struct;
aTmp = cell(length(smoothLTot), length(nCompsTot), length(bwTot), length(smoothKTot), length(nLocConsiderTot), length(regScaleTot));
bestCompDec.dispStore = aTmp;
bestCompDec.dispErrStore = aTmp;
bestCompDec.compChoice = nCompsTot;
bestCompDec.meritFig = aTmp;
bestCompDec.csiErrorTot = aTmp;
bestCompDec.minCsiEstError = aTmp;

bestCompDec.meritFig = aTmp;
bestCompDec.residualEstimationErrorTot = aTmp;
bestCompDec.snrAllSubcarr = aTmp;
bestCompDec.snrCv = aTmp;
bestCompDec.snrAllSub_overall = aTmp;
bestCompDec.snrCv_overall = aTmp;

% trajectory saving
trajOvrSyncTot = cell(length(expOffsetTot));
trajEstShiftTot = cell(length(expOffsetTot));
trajEstSyncTot = cell(length(expOffsetTot));
trajErrSyncTot = cell(length(expOffsetTot));
dispErrSyncTot = cell(length(expOffsetTot));
dispErrSyncNormalizedTot = cell(length(expOffsetTot));

tic
disp_err_mm_tot = [];
traj_err_mm_tot = [];
for iExpOffset = 1:length(expOffsetTot)
for iLayout = 1:length(layoutTot)
for iSnrChoose = 1:length(snrTot)
for iSmoothL = 1:length(smoothLTot)
for iCompChoose = 1:length(nCompsTot)
for iBw = 1:length(bwTot)
for iSmoothK = 1:length(smoothKTot)
for iLocConsider = 1:length(nLocConsiderTot)

%% geenrating data according to input source
% load('tmp');
inputSrc = 'ovr5300'; 
                   % 5300 for intel 5300
                   % 5300_cont for continuous motion experiments with Intel 5300
                   % own for self generated channels
                   % vna for data collected from VNA
                   % vnaKiran for data collected by Kiran
                   % vnaLargeAp for data collected with moving the rx array for each experiment. done only once.
switch inputSrc

case 'ovr5300'
%% CSI traces from Intel 5300
addpath('./cvxgen/') % set to './cvxgenUnix/' if using Unix and './cvxgen/' if using MAC machine
PinPointLocation = '..//cvprAfs/NSDI_Localization/PinPoint123456'; %'/home/mani/PinPoint123456';

expFileOffset = expOffsetTot(iExpOffset);
dispExpDataIdTot = 1:1:4050; %1:2:4850; % 1:8:4900; 
ipIndicesArray = [7 9 10 11];%ipIndices to consider for localization
refIp = 8;

trialIdx = 1;
nTrials = length(trialIdx);

CalibrationFiles = ['log20150107120258.dat'; 'log20150107121905.dat'; 'log20150107124046.dat'; ...
            'log20150107123308.dat'; 'log20150107051313.dat'; 'log20150107122735.dat'; ...
            'log20161009221938.dat'; 'log20160515124858.dat'; 'log20161009221730.dat'; ...
            'log20161009223243.dat'; 'log20161009222529.dat' ];
startFiles = ['log20150107120258.dat'; 'log20150107121905.dat'; 'log20150107124046.dat'; ...
            'log20150107123308.dat'; 'log20150107051313.dat'; 'log20150107122735.dat'; ...
            'log20160415234654.dat'; 'log20160407155621.dat'; ...
            'log20160106223536.dat'; 'log20160407154455.dat'; 'log20160415233954.dat';];
txOrientation = 0; %0;
% APs for which experiments are conducted. this is done for easily copying experimentData here
% CalibrationFiles(ipIndicesArray,:) = [...
% 'log20160515123532.dat';
% 'log20160515124858.dat';
% 'log20160515124025.dat';
% 'log20160515123722.dat';
% ];
% startFiles(ipIndicesArray,:) = [...
% 'log20160525203508.dat';
% 'log20160525204545.dat';
% 'log20160525203431.dat';
% 'log20160525202819.dat';];
% ovr_startFile = '20160525203427.dat';


CalibrationFiles = ['log20150107120258.dat'; 'log20150107121905.dat'; 'log20150107124046.dat'; ...
            'log20150107123308.dat'; 'log20150107051313.dat'; 'log20150107122735.dat'; ...
            'log20161009221938.dat'; 'log20160515124858.dat'; 'log20161009221730.dat'; ...
            'log20161009223243.dat'; 'log20161009222529.dat' ];
startFiles(ipIndicesArray,:) = [...
'log20161011002417.dat';
'log20161011002813.dat';
'log20161011002814.dat';
'log20161011002516.dat';
];
ovr_startFile = '20161011002811.dat';
startFiles(refIp,:) = 'imu20161011002816.log';

%% sub IPS, the subset of IPs you want to consider for tracking
subIps = [1:4];
ipIndicesArray = ipIndicesArray(subIps);

musicComps = [2 2 2 2 ]; % same length as ipIndicesArray

ovrPresent = 1;
if ovrPresent
    % gettingn the OVR orientation data
    ovrDataLoc =  '/Users/mani/Dropbox/ovrTrackdata'; %  '../../../ovrTrackdata'; %
    ovr_filenames = dir( sprintf('%s/*.dat', ovrDataLoc));
    ovr_ExpStartFileIdx = strmatch(ovr_startFile,{ovr_filenames(:).name}) -1 + expFileOffset;
    ovr_filename = ovr_filenames(ovr_ExpStartFileIdx).name
    ovrTmp = load(sprintf('%s/%s', ovrDataLoc, ovr_filename));

    tic
    % reference filename for comparing timestamps of OVR and 5300
    startFile = startFiles(refIp,:);
    filenames = dir( sprintf('%s/ip%d/imu*.log',PinPointLocation,refIp));
    ExpStartFileIdx = strmatch(startFile,{filenames(:).name}) -1 + expFileOffset;
    filename = filenames(ExpStartFileIdx).name;
    refFilename = sprintf('%s/ip%d/capture%s.log', PinPointLocation,refIp,filename(4:end-4));
    toc

    % obtaining orientation of nearest timestamps
    [nearTimeInd, trajOvrSync] = ovrTraj(ovrTmp, refFilename, dispExpDataIdTot, 0);
    [yaw, pitch, roll] =  ovrOrient(ovrTmp, nearTimeInd);

    % find the indices with orientation more or less equal to the start orientation
    smallOriThres = 0.005;
    orientConsiderInd = find( abs(diff(yaw))>smallOriThres, 1, 'first' ) : find( abs(diff(yaw))>smallOriThres, 1, 'last' );
    orientConsider = yaw(orientConsiderInd);
    orientThres = 30; %9;
    histEdges = linspace(-180, 180, 1 + 360/(2*orientThres));
    [binCount,histEdges] = histcounts(orientConsider, histEdges);
    [~, maxBin] = max(binCount);
    sameOrientIndTmp = find(yaw>histEdges(maxBin) & yaw<histEdges(maxBin+1));
    sameOrientInd = dispExpDataIdTot((sameOrientIndTmp));
    sprintf('number of points with more or less same orientation is %d', length(sameOrientInd))
    ExpDataIdTot = sameOrientInd(round( linspace(1,length(sameOrientInd),min(600,length(sameOrientInd))) )); %sameOrientInd(min(1, end-500):end);
    % 
else
    ExpDataIdTot = dispExpDataIdTot;
    yaw = zeros(size(dispExpDataIdTot));
end

nIpsTot = length(ipIndicesArray); % total number of APs used for the experiment. If 11 APs are there in total, then nIpsTot = 11;

if ~exist('wiredCalib', 'var')
    wiredCalib = RxResponse3Tx3Rx(CalibrationFiles, ipIndicesArray);
end

csi_trace_all_multTrials = cell(1,nIpsTot,length(ExpDataIdTot));
csi_trace_disp_all_multTrials = cell(1,nIpsTot,length(dispExpDataIdTot));
for ipIndex = 1:length(ipIndicesArray)
    ipNo = ipIndicesArray(ipIndex);

    startFile = startFiles(ipNo,:);
    filenames = dir( sprintf('%s/ip%d/log*.dat',PinPointLocation,ipNo));
    ExpStartFileIdx = strmatch(startFile,{filenames(:).name}) -1 + expFileOffset;
    filename = filenames(ExpStartFileIdx).name
    if ipNo == refIp
        refFilename = filename;
    end
    csi_trace = read_bf_file(sprintf('%s/ip%d/%s',PinPointLocation,ipNo,filename)); % replace with location tmp.dat
    sprintf(' has %d packets', length(csi_trace))

    csi_trace_all_multTrials = cellFromCsi(csi_trace_all_multTrials, csi_trace, wiredCalib, ipNo, ipIndex, ExpDataIdTot, M, N, T_3Tx, SubCarrInd);
    csi_trace_disp_all_multTrials = cellFromCsi(csi_trace_disp_all_multTrials, csi_trace, wiredCalib, ipNo, ipIndex, dispExpDataIdTot, M, N, T_3Tx, SubCarrInd);
end

nRx = nIpsTot;
nSegments = length(ExpDataIdTot);
nSegmentsDisp = length(dispExpDataIdTot);

if length(nLocConsiderTot)==1
    nLocConsiderTot = nSegments;
end


end



%% find the AoD 
iTrial = 1;
csi_trace_all = squeeze(csi_trace_all_multTrials(iTrial,:,:));
csi_trace_disp_all = squeeze(csi_trace_disp_all_multTrials(iTrial,:,:));
if nRx == 1
    csi_trace_all = squeeze(csi_trace_all_multTrials(iTrial,:,:)).';
    csi_trace_disp_all = squeeze(csi_trace_disp_all_multTrials(iTrial,:,:)).';
end
clearvars csi_trace_all_multTrials csi_trace_disp_all_multTrials
% csi_trace_all_Full = squeeze(csi_trace_all_Full_multTrials(1,:,:)).';

aoaDodAllSegments = cell(nRx,1);
% paramStore_3d_2pkts = cell(nRx, nSegmentsValid);
% dispEstFromDodRayStore = cell(nRx, nSegmentsValid);
% dispEstMultPktStore = cell(nRx, nSegmentsValid);
% % initialization for the AoA-DoD case
aoaDodOnly = cell(nRx, N*nSegmentsDisp);
% csi_aoaDodOnly = cell(nRx, N*nSegmentsDisp);
for iRx = 1:nRx
    nExpsForSegments = floor(nSegments/nLocConsiderTot(iLocConsider));
    for expEndLoc = 1:nExpsForSegments
        nApTogether = 1;
        nTxTogether = nLocConsiderTot(iLocConsider); % nSegments before
        txGap = 1;
        expConsider = (expEndLoc-1)*nLocConsiderTot(iLocConsider)+1;

      for iExp = expConsider
        paramRange = struct;
        paramRange.delayRange = [0 170]*1e-9;
        paramRange.angleRange = 90*[-1 1];
        paramRange.GridPts = [1 101 101];
        paramRange.K = M; %smoothKTot(iSmoothK);  
        paramRange.L = 1; %floor(N*smoothLTot(iSmoothL));
        paramRange.T = T_3Tx;

        paramRange.circularTx = 1;
        paramRange.dTx = dTx;
        txMatParams = struct;
        txMatParams.dTx = dTx;
        txMatParams.circular = 1;
        paramRange.deltaRange = [0 359]; %(dTx)*[-1 1];
        paramRange.generateAtot = 1;

        vna_response = [];
        for iTx = 1:nTxTogether
            vna_response_tmp = formatCSI(csi_trace_all{iRx, iExp+(iTx-1)*txGap}, N, M, T_3Tx, paramRange.K, paramRange.L, paramRange.T);
            vna_response = [vna_response vna_response_tmp];
        end
        paramRange.X = vna_response;
        maxRapIters = Inf;
        useNoise = 0;
        do_second_iter = 0;
        paramRange.seconditerGridPts = [1 51 21];
        paramRange
        trueParams.txOrientation = txOrientation;
        trueParams.nComps = min([musicComps(iRx) size(paramRange.X)-1]);    
        trueParams.nCompsDisp = trueParams.nComps;    

        if ~isfield(trueParams,'nComps') 
            aoaDodAllSegmentsTmp = backscatterEstimationMusic(ones(N*M*T_3Tx,1), M, N, c, fc,...
            T_3Tx, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter);
        else
            aoaDodAllSegmentsTmp = backscatterEstimationMusic(ones(N*M*T_3Tx,1), M, N, c, fc,...
            T_3Tx, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(trueParams.nComps))   
        end

        for iTx = 1:nTxTogether
            aoaDodAllSegments{iRx,iExp + (iTx-1)*txGap} = aoaDodAllSegmentsTmp;
        end

      end
    end

end
clearvars csi_trace_all
%% finding the phase of the first tx-rx antenna pair at each of the lcations
tic
ovrYawDir = 0; -1;
expStruct = v2struct(fc,M,fs,c,d,dTx,SubCarrInd,N,fgap,T_3Tx,tofInd,aoaInd,dodInd,alphaInd);
for iRx = 1:nRx
    for iExp = 1:size(csi_trace_disp_all,2)
        csi_trace_exp = reshape(csi_trace_disp_all{iRx,iExp}, N, M, T_3Tx);
        for iSubCarr = 1:N
            aoaDodAllSegmentsTmp = aoaDodAllSegments{iRx, 1};
            aoaDodAllSegmentsTmp(:,dodInd) = aoaDodAllSegmentsTmp(:,dodInd) - ovrYawDir*(yaw(iExp)-yaw(1)); % modifying the the DoD to correct for yaw
            Ahat = steeringMatBackscatter(fc, T_3Tx, aoaDodAllSegmentsTmp(:,dodInd), M, aoaDodAllSegmentsTmp(:, aoaInd), c, ...
                    [0], fgap, aoaDodAllSegmentsTmp(:, tofInd)*1e-9, d, txMatParams); 
            csiVecTmp = reshape(csi_trace_exp(iSubCarr,:,:), [],1);
            aoaDodAllSegmentsTmp(:,alphaInd) = Ahat\csiVecTmp;
            aoaDodOnly{iRx, (iExp-1)*N + iSubCarr} = aoaDodAllSegmentsTmp; 
            % csi_aoaDodOnly{iRx, (iExp-1)*N + iSubCarr} = csiVecTmp;
            % initResid((iExp-1)*N + iSubCarr) = norm(csiVecTmp - Ahat*aoaDodAllSegmentsTmp(:,alphaInd))^2/norm(csiVecTmp)^2;
        end
    end
    [aoaDodOnly{iRx,1}(:,1:3) db(abs(aoaDodOnly{iRx,1}(:,4)))] 
end

toc
clearvars csi_trace_disp_all
%% finding the displacements
[trajEst, dispEst_fromalpha] = dispOvr5300(aoaDodOnly,expStruct, nSegmentsDisp, nRx, txOrientation);
figure(5);
plot( trajEst(:,1), trajEst(:,2), 'd-', 'LineWidth', 2); axis equal; 
xlabel('mm'); ylabel('mm')
sum(norms(diff(trajEst),2,2))

% f4 = figure(6);
% for iPt = 1:10:size(trajEst,1)
%     plot( trajEst(1:iPt,1), trajEst(1:iPt,2), 'd-', 'LineWidth', 2); axis equal; 
%     xlabel('mm'); ylabel('mm')
%     waitforbuttonpress
%     %pause(0.1)
% end

% trajEstRotate = trajEst;
% trajEstRotate(:,2) = -trajEstRotate(:,2);

trajEstSync = trajEst;
[d,trajEstShift] = procrustes(trajOvrSync, trajEstSync, 'scaling', false, 'reflection', false);
%%
f2 = figure(5);
plot( trajEstShift(:,2), -trajEstShift(:,1), '-', 'LineWidth', 2); hold on; 
plot( trajOvrSync(:,2), -trajOvrSync(:,1), ':', 'LineWidth', 4);  hold off; axis equal; axis tight; 
legend('WiCapture', 'Oculus Rift', 'Location','NorthOutside','Orientation','horizontal')
xlabel('motion along x-axis in mm'); ylabel('motion along y-axis in mm')
% set(gca,'position',[0 0 1 1],'units','normalized')
% saveas(f2,sprintf('../../PinPoint123456/ovrWiFiFigs/trajOvrWifi%d.png',expFileOffset))
% saveas(f2,sprintf('/Users/mani/Dropbox/ovrTraces/trajOvrWifi%d.png',expFileOffset))

% %% tikz picture for fellowship essay
% addpath_recurse('./export_fig_fold')
% addpath_recurse('./matlab2tikzFold')
% tikzDest = '/Users/mani/Documents/fbFellowshipApplication/msrAppLatex';
% cleanfigure; 
% % matlab2tikz(  [tikzDest '/' 'trace.tikz'],'interpretTickLabelsAsTex',true, 'width','0.8\textwidth', 'figurehandle',f2);
% matlab2tikz( [tikzDest '/' 'trace.tikz'],'interpretTickLabelsAsTex',true, 'width','0.8\textwidth', 'extraAxisOptions',['ylabel shift={-7pt},', 'xlabel shift={-2pt},'...
%     'legend style={font=\large},', 'ylabel style={font=\large},', 'xlabel style={font=\large},', 'xticklabel style={font=\small},', 'yticklabel style={font=\small}'])

%%
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

trajOvrSyncTot{expFileOffset} = trajOvrSync;
trajEstShiftTot{expFileOffset} = trajEstShift;
trajEstSyncTot{expFileOffset} = trajEstSync;
trajErrSyncTot{expFileOffset} = trajErrSync;
dispErrSyncTot{expFileOffset} = dispErrSync;
dispErrSyncNormalizedTot{expFileOffset} = dispErrSyncNormalized;

% f3 = figure(7);
% for iPt = 1:10:size(trajEstShift,1)
%     plot( trajEstShift(1:iPt,1), trajEstShift(1:iPt,2), 'd-', 'LineWidth', 2); hold on; 
%     plot( trajOvrSync(1:iPt,1), trajOvrSync(1:iPt,2), '-', 'LineWidth', 4);  hold off; axis equal; axis tight; 
%     legend('from WiFi', 'true path', 'Location','NorthOutside','Orientation','horizontal')
%     xlabel('mm'); ylabel('mm')
%     waitforbuttonpress
%     %pause(0.1)
% end

end
end
end
end
end
end
end
% clearvars  csi_trace_all_multTrials csi_trace csi_trace_all csi_aoaDodOnly  paramRange vna_response paramEachSubcarr filenames ovr_filenames aoaDodOnly
% % save(sprintf('../../PinPoint123456/ovrWiFiFigs/trajOvrWiFi%d.mat',expFileOffset))
% save(sprintf('/Users/mani/Dropbox/ovrTraces/trajOvrWiFi%d.mat',expFileOffset))
% save(sprintf('../../PinPoint123456/ovrWiFiFigs/trajOvrWiFi%d.mat',expFileOffset))
end

