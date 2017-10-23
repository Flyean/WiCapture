function trajEst = wicaptureAlgo(csi_trace_all, csi_trace_disp_all, segStart, segEnd,nSegmentsDisp,ExpDataIdTot,dispExpDataIdTot,yaw,yaw_ExpDataId)
% % % INPUT:
% csi_trace_all: (number of receivers x number of WiFi packets) dimensional cell which contains CSI of packets. These packets are used for estimating propagation path parameters
% csi_trace_disp_all: (number of receivers x number of WiFi packets) dimensional cell which contains CSI of packets. These packets are for which position of the WiFi transmitter is estimated. 
% segStart: multiple packets are transmitted. The packets are grouped into length(segStart) number of groups. segStart is the starting packet index for each of the groups.
% segEnd: segEnd is the ending packet index for each of the groups.
% nSegmentsDisp: number of packets for which position needs to be estimated.
% ExpDataIdTot: A large number of packets are transmitted for each experiment. CSI in csi_trace_all are CSI from pakets indexed by ExpDataIdTot.
% dispExpDataIdTot: A large number of packets are transmitted for each experiment. CSI in csi_trace_disp_all are CSI from pakets indexed by dispExpDataIdTot.
% yaw: yaw of the WiFi transmitter when the transmitter is transmitting packets index by dispExpDataIdTot
% yaw_ExpDataId: yaw of the WiFi transmitter when the transmitter is transmitting packets index by ExpDataIdTot

fc = 5.63e9; % center frequency
M = 3;  % number of antennas at the receiver
fs = 40e6; % WiFi channel bandwidth
c = 3e8; % speed of light
d = 2.6e-2; % distance between adjacent antennas at the receiver
dTx = 2.6e-2; % distance between adjacent antennas at the transmitter
SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % subcarriers at which CSI is available
N = length(SubCarrInd); % number of subcarriers where CSI is available
fgap = 312.5e3; % frequency gap between adjacent WiFi subcarriers
T_3Tx = 3; % number of transmitter antennas
% The propagation path parameters are stored as a matrix where different columns correspond to different parameters
tofInd = 1; % column index which corresponds to ToF (time of flight)/elevation angle of the paths
elevDodInd = 1; % column index which corresponds to ToF/elevation angle of the paths
aoaInd = 2; % column index which corresponds to AoA (angle of arrival) of the paths
dodInd = 3; % column index which corresponds to DoD (direction of departure) of the paths
alphaInd = 4; % column index which corresponds to complex attenuation of the paths
ovrYawDir = -1;

nRx = 4;  % number of receivers/APs (access points) at which CSI is collected
musicComps = [2 2 2 2]; % number of paths between transmitter and each of the four receivers


nContSegments = length(segStart); % multiple packets are transmitted. The packets are 
% are grouped into length(segStart) number of groups. The first group contains segStart(1)
% to segStart(end) packets. nContSegments is the number of groups the packets from the experiment are divided into.

    aoaDodAllSegments = cell(nRx,1); % cell structure storing the estimated propagation path parameters
    aoaDodOnly = cell(nRx, N*nSegmentsDisp); % cell structure storing the estimated propagation path parameters for each packet at each AP
    csi_aoaDodOnly = cell(nRx, N*nSegmentsDisp);
    meanYaw = Inf*ones(nContSegments,1); % mean of the yaw for all the packets in a particular group
    
%%    % % % % % % find the propagation path parameters
    for iRx = 1:nRx
        for iContSegment = 1:nContSegments
          %  find the packets for which orientation of the transmitter WiFi chip does not change much.         
          orientThres = 10;
          [sameOrientTxInd, meanYaw(iContSegment)] = findSameOrient(yaw_ExpDataId(segStart(iContSegment):segEnd(iContSegment)),orientThres);
          
          txGap = 1;
          % find the proapgation path parameters for each group of packets
          for iExp = segStart(iContSegment)
            paramRange = struct; % structure to store parameters required for applying MUSIc algorithm
            paramRange.delayRange = [0 170]*1e-9; % range of ToF values considered in MUSIC 
            paramRange.angleRange = 90*[-1 1]; % range of AoA values considered in MUSIC 
            paramRange.elevationRange = [90 90]; % range of AoD (angle fo departure elevaion) values considered in MUSIC 
            paramRange.GridPts = [1 101 101 1+floor(diff(paramRange.elevationRange))]; % number of grid points for [ToF AoA AoD AoD_elevation] parameters
            do_second_iter = 0;
            paramRange.seconditerGridPts = [1 51 21 21];
            paramRange.K = M; % number of receive antennas in 'smoothed' CSi matrix
            paramRange.L = 1; % number of subcarriers in 'smoothed' CSi matrix
            paramRange.T = T_3Tx; % number of transmit antennas in 'smoothed' CSi matrix

            paramRange.circularTx = 1; % 1 if the transmit antennas form a circular array
            paramRange.dTx = dTx;  % distance between transmit antennas
            txMatParams = struct; 
            txMatParams.dTx = dTx;  % distance between transmit antennas
            txMatParams.circular = 1; % 1 if the transmit antennas form a circular array
            paramRange.deltaRange = [0 359]; % range of AoD (angle of departure) values considered in MUSIC 
            paramRange.generateAtot = 0; % 0 if you want to use pre-computed and stored steering vector values
            
            % % Form X by concatenating matrices from different packets
            % % Equation 7 in the paper            
            vna_response = [];
            for iTxTmp = 1:length(sameOrientTxInd)
                iTx = sameOrientTxInd(iTxTmp);
                vna_response_tmp = formatCSI(csi_trace_all{iRx, iExp+(iTx-1)*txGap}, N, M, T_3Tx, paramRange.K, paramRange.L, paramRange.T);
                vna_response = [vna_response vna_response_tmp];
            end
            
            paramRange.X = vna_response; % input matrix to the MUSIC algorithm 
            maxRapIters = Inf;
            useNoise = 0;
            trueParams.nComps = min([musicComps(iRx) size(paramRange.X)-1]);    
            trueParams.nCompsDisp = trueParams.nComps;    

            if iContSegment == 1
                % apply MUSIC algorithm
                    aoaDodAllSegmentsTmp = backscatterEstimationMusic3D(ones(N*M*T_3Tx,1), M, N, c, fc,...
                    T_3Tx, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(trueParams.nComps))   
                
                % % % aoa dod obtained by using all the packets
                aoaDodAllCsi = aoaDodAllSegmentsTmp;
            end
            
            % % % trick to reduce the computation needed for groups of packets other than the first group
            % % % instead of sweeping through entire angle range, just sweep through +-10 degrees of AoA, AoD estimate of previous group
            if iContSegment > 1        
                paramRange.X = vna_response; % resetting the CSI responses to be considered just to the small segment instead of all the responses
                aoaDodAllSegmentsTmp = aoaDodAllSegments{iRx,iExp-1};
                for iMp = 1:size(aoaDodAllSegments{iRx,iExp-1},1)
                        % % if you want to consider segment by segment
                        paramRange.angleRange = aoaDodAllSegments{iRx,iExp-1}(iMp, aoaInd) + 10*[-1 1];
                        paramRange.elevationRange = aoaDodAllSegments{iRx,iExp-1}(iMp, elevDodInd) + 0*[-1 1];
                        paramRange.deltaRange = aoaDodAllSegments{iRx,iExp-1}(iMp, dodInd) + 10*[-1 1]...
                                                    - ovrYawDir*(meanYaw(iContSegment)-meanYaw(iContSegment-1)) ;
                        paramRange.GridPts = [1 21 21 1]; %[1 51 51 51];
                    aoaDodSmallSegment = backscatterEstimationMusic3D(ones(N*M*T_3Tx,1), M, N, c, fc,...
                        T_3Tx, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, ones(trueParams.nComps));
                    aoaDodAllSegmentsTmp(iMp,:) = aoaDodSmallSegment(1,:);
                end
            end

            % saving the propagation path parameters like AoA AoD
            for iTx = 1:(segEnd(iContSegment) - segStart(iContSegment) + 1)
                aoaDodAllSegments{iRx,iExp + (iTx-1)*txGap} = aoaDodAllSegmentsTmp;
            end

          end
        end

    end
    
    %% finding the complex attenuation along each of the paths
    tic
    % ExpDataIdTot is the indices of the WiFi packets using which WiCapture calculates the propagation path parameters
    % dispExpDataIdTot is the indices of the WiFi packets at which WiCapture wants to find the position of the WiFi transmitter
    smallestExpEst = knnsearch(ExpDataIdTot(:), dispExpDataIdTot(:)); 
    
    expStruct = v2struct(fc,M,fs,c,d,dTx,SubCarrInd,N,fgap,T_3Tx,tofInd,aoaInd,dodInd,alphaInd,elevDodInd);
    for iRx = 1:nRx

        for iExp = 1:size(csi_trace_disp_all,2)
            csi_trace_exp = reshape(csi_trace_disp_all{iRx,iExp}, N, M, T_3Tx);
            
            aoaDodAllSegmentsTmp = aoaDodAllSegments{iRx, smallestExpEst(iExp)};

            do2DCorr = 1;
            if do2DCorr
               aoaDodAllSegmentsTmp(:,dodInd) = aoaDodAllSegmentsTmp(:,dodInd) - ovrYawDir*(yaw(iExp)-meanYaw(iContSegment)); % modifying the the DoD to correct for yaw 
            end

                Ahat = [];
                for iComp = 1:size(aoaDodAllSegmentsTmp,1)
                    u_s = (d*fc/c)*sin(deg2rad(aoaDodAllSegmentsTmp(iComp,aoaInd)));
                    Ahat(:,iComp) = circularGridSampleBackscatter3D(fc, [], [], M, u_s, c, 0, fgap, 0, dTx,...
                        aoaDodAllSegmentsTmp(iComp,dodInd), aoaDodAllSegmentsTmp(iComp,elevDodInd) );
                end
                
            for iSubCarr = 1:N
                csiVecTmp = reshape(csi_trace_exp(iSubCarr,:,:), [],1);
                aoaDodAllSegmentsTmp(:,alphaInd) = Ahat\csiVecTmp;
                aoaDodOnly{iRx, (iExp-1)*N + iSubCarr} = aoaDodAllSegmentsTmp; 
                csi_aoaDodOnly{iRx, (iExp-1)*N + iSubCarr} = csiVecTmp;
            end
        end
        % displaying the propagation path parameters        
        [aoaDodOnly{iRx,1}(:,1:3) db(abs(aoaDodOnly{iRx,1}(:,4)))] 
    end
    toc
    tic
    aoaDodOnlyOrig = aoaDodOnly;
    toc
    
   
%% findidng the displacement of the WiFi transmitter between successive packets in dispExpDataIdTot indices
txOrientation = 0; % orientation of the transmiter for different packets
tic
[trajEst, dispEst_fromalpha, A_dod_comb_allRxSave, compMinVar_dispDiffMusic_comb_allRxSave] = dispOvr5300(aoaDodOnly,expStruct, nSegmentsDisp, nRx, txOrientation);
% [trajEst, dispEst_fromalpha] = dispOvr5300_3D(aoaDodOnly,expStruct, nSegmentsDisp, nRx, txOrientation,elevDodInd);
toc

clearvars  csi_trace_all_multTrials csi_trace csi_trace_all csi_aoaDodOnly  paramRange vna_response paramEachSubcarr aoaDodOnlyOrig


