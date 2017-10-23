%% estimates the propagation path parameters of all the paths from CSI
function [Parameters, varargout] = backscatterEstimationMusic3D(csi_trace, M, N, c, fc, Ttot, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, varargin)
% INPUTS:
% csi_trace : (N*M*Ttot)x1 vector containing the csi trace
% M: number of receive antennas
% Ttot: number of transmit antennas
% fc: center frequency
% c: speed of light
% do_second_iter: use 0
% delayRange: 1x2 dimensional vector with the minimum and maximum values of the ToF in s for which MUSIC spectrum needs to be calculated.
% SubCarrInd: indices of subcarriers where CSI is provided
% N: number of such subcarriers
% fgap: frequency gap between successive subcarriers in Hz
% d: distance between adjacent antennas in m
% paramRange: structure containing parameters for calculating MUSIC spectrum. See wicaptureAlgo.m for details.
% maxRapIters: maximum number of paths that you want to consider. Use Inf
% useNoise: use 0
% do_second_iter: use 0;

% OUTPUTS:
% Parameters : Each row of parameters contains [ToF(ns)/Elevation_angle(degrees) AoA(degrees) displacement(mm)/azimuth_angle(degrees)]


% one of the heuristics to decide the number of paths
nComps = [];
corrThr = 0.4; % correlation threshold to decide after how many components we decide to cutoff the MUSIC process
                % previously used a value of 0.83 . This value should be square of cosd of cutoff angle with the projection plane      
if ~isempty(varargin)
    parametersTrue = varargin{1};
    nComps = size(parametersTrue,1);
    if nComps==0, nComps = []; end;
    corrThr = 0;
end
%% MUSIC parameters
% % % For AoA-ToF smoothing
% K = floor(M/2)+1; % K is the number fo subset of antennas chosen for smoothing music. Default value is floor(M/2)+1
% L = floor(N/2); % L is the number fo subset of subcarriers chosen for smoothing music. Default value is floor(N/2)
% % % For mobile localization 
if ~isfield(paramRange,'K')
    K = floor(M/2)+1; 
    L = floor(N/2);
    T = floor(Ttot/2)+1;
else
    K = paramRange.K; 
    L = paramRange.L;
    T = paramRange.T;
end

% do_second_iter = 0; % do two iterative gridding in MUSIC
% deltaRange: 1x2 dimensional vector with the minimum and maximum values of the AoD azimuth angle in degrees for which MUSIC spectrum needs to be calculated.
% delayRange: 1x2 dimensional vector with the minimum and maximum values of the ToF in s for which MUSIC spectrum needs to be calculated.
% elevationRange: 1x2 dimensional vector with the minimum and maximum values of the AoD elevation angle in degrees for which MUSIC spectrum needs to be calculated.
% angleRange: same as above but for AoA in degrees
% GridPts: 1x3 vector containing number of grid points for [ToF AoA AoD] parameters
if ~isfield(paramRange,'delayRange')
    % % % % FIXED DELAY, DELTA RANGES USED IN PREVIOUS VERSION
    delayRange = [-25 25]*1e-9; % [-25 70]
    deltaRange = [-c/2/fc  c/2/fc];
    angleRange = [-90 90];
    GridPts = [100 100 100];
else
    % % % % TAKING DELAY, DELTA, ANGLE ranges from INPUT
    delayRange = paramRange.delayRange;
    deltaRange = paramRange.deltaRange;
    angleRange = paramRange.angleRange;
    GridPts = paramRange.GridPts;
end
MaxAngle = angleRange(2);
MinAngle = angleRange(1);
if isempty(maxRapIters)
    maxRapIters = Inf;
end

% setting the grid points
if do_second_iter
    if ~isfield(paramRange,'seconditerGridPts') %if isempty(paramRange)
        GridPts = [70 70 35];  % [100 50 50] previously used % GridPts correspond to ToF, AoA and delta respectively
        MaxAngle = MaxAngle*(GridPts(2)-1)/(GridPts(2)+1);
        MinAngle = -MaxAngle;
        seconditerGridPts = [15 25 2]; %[15 25 5];
    else
        seconditerGridPts = paramRange.seconditerGridPts;
    end
else
%     GridPts = paramRange.GridPts; % [300 5 1]; % GridPts correspond to ToF, AoA and delta respectively %[1 180 90]; %
%     MinAngle = -90;
%     MaxAngle = 90;
    seconditerGridPts = [];
end


%% parameters for circular array
if ~isfield(paramRange,'circularTx')
    paramRange.circularTx = 0;
end


if paramRange.circularTx == 1
    
    if ~isfield(paramRange, 'deltaRange')
       deltaRange = 0:359;
    else
       deltaRange = paramRange.deltaRange;
    end
    
    if ~isfield(paramRange, 'elevationRange')
       elevationRange = 0:90;
    else
       elevationRange = paramRange.elevationRange;
    end

    dTx = paramRange.dTx;
   
end


%%     Vector containing the components of each ray at each grid point
[aTot,GridStart,GridSpacing, delayGridValue, u_sGridValue, deltaGridValue, elevationGridValue] = gridVecBackscatter3D(deltaRange, paramRange.K, T, d, fc, c, do_second_iter, delayRange, SubCarrInd, N, fgap, GridPts, MaxAngle, MinAngle, elevationRange, paramRange.generateAtot);
%     aTot = gpuArray(aTot);

EigDiffCutoff = 4; % previously used value 4 % difference between successive eigen values to consider the next eigenvalue as noise


%% formatting the input CSI matrix to provide input to MUSIC algorithm
if ~isfield(paramRange,'X')
    X = formatCSI(csi_trace, N, M, Ttot, K, L, T);
else
    X = paramRange.X;
end


%% Applying MUSIC 
[Pn,Ps,Qn,Qs,EigenInfo] = GetQnBackscatter(X,EigDiffCutoff, nComps);

delayFromMusic = [];
angleFromMusic = [];
deltaFromMusic = [];
elevationFromMusic = [];
if ~useNoise
    nIters = min(maxRapIters,size(Qs,2)); % number of iterations for RAP MUSIC
else
    nIters = 1; 
end
maxCorr = zeros(nIters,1);
doBreak = 1; % whether to apply break statement
      
if paramRange.circularTx == 0  % if the transmitter has a uniform linear array

        for k=1:nIters
            % DelayStartStop = DelayStartStopOrig;%+(k-1)*2e-9; 
            % [delayFromMusicTmp,angleFromMusicTmp, maxCorr(k)]= SingleSnapshotRAPMusic(Qn,Qs,fc,fgap,d,K,L,delayFromMusic,angleFromMusic,SubCarrInd,DelayStartStop);
            [delayFromMusicTmp,angleFromMusicTmp, deltaFromMusicTmp, maxCorr(k),music_spectrum] = RAPMusicGridMaxBackscatter(aTot,GridStart,GridSpacing,GridPts,Qn,Qs,fc,fgap,d,K,L,delayFromMusic,angleFromMusic, deltaFromMusic, ...
                                        SubCarrInd, deltaGridValue, u_sGridValue, delayGridValue, T, c, do_second_iter, seconditerGridPts, useNoise);
            if ~doBreak
                delayFromMusic = [delayFromMusic; delayFromMusicTmp];
                angleFromMusic = [angleFromMusic; angleFromMusicTmp];
                deltaFromMusic = [deltaFromMusic; deltaFromMusicTmp];
            else
                % % % % % % %             Previous way of doing RAP MUSIC where I have a
                % % % break statement to break when the maximum correlation goes below certian
                % % % value
                if k==1
                    delayFromMusic = [delayFromMusic; delayFromMusicTmp];
                    angleFromMusic = [angleFromMusic; angleFromMusicTmp];
                    deltaFromMusic = [deltaFromMusic; deltaFromMusicTmp];
                else
                    if maxCorr(k)>corrThr*max(maxCorr)
                        delayFromMusic = [delayFromMusic; delayFromMusicTmp];
                        angleFromMusic = [angleFromMusic; angleFromMusicTmp]; 
                        deltaFromMusic = [deltaFromMusic; deltaFromMusicTmp];
                    else
                        break % In nested loops, break exits only from the for/while loop in which it occurs.
                    end
                end
            end
        end

        if ~doBreak
            allParameters = [delayFromMusic*1e9 angleFromMusic deltaFromMusic*1e3];
            varargout{1} = allParameters;
            % if you do not break, then consider only those parameters which give
            % significant MUSIC spectrum value
            cutoffEntry = find(maxCorr<corrThr*maxCorr(1),1)-1; % Find the first entry that is less than 0.83 times the first maxCorr entry
            if isempty(cutoffEntry)
               cutoffEntry = length(maxCorr); 
            end
            delayFromMusic = delayFromMusic(1:cutoffEntry);
            angleFromMusic = angleFromMusic(1:cutoffEntry);
            deltaFromMusic = deltaFromMusic(1:cutoffEntry);
        end

end


%% for ciurcular array, finding the parameers of interest
if paramRange.circularTx == 1 % if the transmitter has a circular array of 3 antennas

        for k=1:nIters
            % DelayStartStop = DelayStartStopOrig;%+(k-1)*2e-9; 
            % [delayFromMusicTmp,angleFromMusicTmp, maxCorr(k)]= SingleSnapshotRAPMusic(Qn,Qs,fc,fgap,d,K,L,delayFromMusic,angleFromMusic,SubCarrInd,DelayStartStop);
            [delayFromMusicTmp,angleFromMusicTmp, deltaFromMusicTmp, elevationFromMusicTmp, maxCorrTmp,music_spectrum] = circularRAPMusicGridMaxBackscatter3D(aTot,GridStart,GridSpacing,GridPts,Qn,Qs,fc,fgap,d,K,L,delayFromMusic,angleFromMusic, deltaFromMusic, elevationFromMusic, ...
                                        SubCarrInd, deltaGridValue, u_sGridValue, delayGridValue, elevationGridValue, T, c, do_second_iter, seconditerGridPts, useNoise, dTx);
            if ~doBreak
                delayFromMusic = [delayFromMusic; delayFromMusicTmp];
                angleFromMusic = [angleFromMusic; angleFromMusicTmp];
                deltaFromMusic = [deltaFromMusic; deltaFromMusicTmp];
                elevationFromMusic = [elevationFromMusic; elevationFromMusicTmp];
            else
                % % % % % % %             Previous way of doing RAP MUSIC where I have a
                % % % break statement to break when the maximum correlation goes below certian
                % % % value
                if k==1
                    delayFromMusic = [delayFromMusic; delayFromMusicTmp];
                    angleFromMusic = [angleFromMusic; angleFromMusicTmp];
                    deltaFromMusic = [deltaFromMusic; deltaFromMusicTmp];
                    elevationFromMusic = [elevationFromMusic; elevationFromMusicTmp];
                else
                    maxCorr(k) = maxCorrTmp;
                    if maxCorr(k)>corrThr*max(maxCorr)
                        delayFromMusic = [delayFromMusic; delayFromMusicTmp];
                        angleFromMusic = [angleFromMusic; angleFromMusicTmp]; 
                        deltaFromMusic = [deltaFromMusic; deltaFromMusicTmp];
                        elevationFromMusic = [elevationFromMusic; elevationFromMusicTmp];
                    else
                        break % In nested loops, break exits only from the for/while loop in which it occurs.
                    end
                end
            end
        end
        deltaFromMusic = deltaFromMusic*1e-3; % divde by 1000 as it gets multiplied by 1000 later.
end

%% finding the amplitude of each component
alphaFromMusic = zeros(length(delayFromMusic),1);
if paramRange.circularTx == 1
    if ~isfield(paramRange, 'X')
        Ahat = [];
        for compNo = 1:length(delayFromMusic)
            u_s = (d*fc/c)*sind(angleFromMusic(compNo));
            Ahat(:,compNo) = circularGridSampleBackscatter3D(fc, Ttot, deltaFromMusic(compNo), M, u_s, c, SubCarrInd(1:N), fgap, delayFromMusic(compNo), dTx, deltaFromMusicPrev(compNo), elevationFromMusicPrev(compNo) );
        end
        alphaFromMusic = Ahat\csi_trace;
        residualEstimationError = norm(csi_trace - Ahat*alphaFromMusic)/norm(csi_trace);
        varargout{1}.residualEstimationError = residualEstimationError;
        % alphaRealImg = [real(Ahat) -imag(Ahat); imag(Ahat) real(Ahat)]\[real(csi_trace); imag(csi_trace)];
        % alphaRealImg2Col = reshape(alphaRealImg, [], 2);
        % alphaComplex = alphaRealImg2Col(:,1) + 1i*alphaRealImg2Col(:,2)
        % residualEstimationError = norm(csi_trace - Ahat*alphaComplex)/norm(csi_trace)
    end
end
% [~,alphFromMusic,ResidualNorm] = GetResidual(delayFromMusic,angleFromMusic,vna_response,d,lambda,fgap,SubCarrInd);
% ParametersDisplay = [delayFromMusic*1e9 angleFromMusic db(abs(alphFromMusic))]
% % % % % % % % %      ns           degrees             mm           normal
Parameters = [elevationFromMusic angleFromMusic deltaFromMusic*1e3 abs(alphaFromMusic) angle(alphaFromMusic)];

% chirp_z = 1;
% if chirp_z == 1
%     Parameters(:,2) = asind(Parameters(:,2));
% end