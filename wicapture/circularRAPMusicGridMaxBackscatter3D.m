%% calculates the MUSIC spectrum and obtains the propagation path parameters of different paths
function [delayFromMusic, angleFromMusic, deltaFromMusic, elevationFromMusic, maxCorr, music_spectrum_plot] = circularRAPMusicGridMaxBackscatter3D(aTot,GridStart,GridSpacing,GridPts,Qn,Qs,fc,fgap,d,K,L,delayFromMusicPrev,angleFromMusicPrev, deltaFromMusicPrev, elevationFromMusicPrev,...
                                                    SubCarrInd, deltaGridValue, u_sGridValue, delayGridValue, elevationGridValue, T, c, do_second_iter, seconditerGridPts, useNoise, dTx)
% % % INPUT:
% aTot: matrix with each column being a steering vector
% GridStart: 1x4 vector: smallest value of each of [ToF AoA AoD AoD_elevation] parameters at which MUSIC spectrum is estimated.
% GridSpacing: 1x4 vector: difference between successive grid values of [ToF AoA AoD AoD_elevation] parameters at which MUSIC spectrum is estimated.
% GridPts: 1x4 vector: number of points along each dimensions of [ToF AoA AoD AoD_elevation] parameters at which MUSIC spectrum is estimated.
% Qn: matrix where each column is a basis vector of noise subspace (read MUSIC algorithm)
% Qs: matrix where each column is a basis vector of signal subspace (read MUSIC algorithm)                                               
% fc: center frequency
% fgap: frequency gap between subcarriers
% d: distance between receiver antennasin m
% K: number of receive antennas in 'smoothed CSI'
% L: number of subcarriers in 'smoothed CSI'
% delayFromMusicPrev: (number of paths calculated so far)-dimensional vector containing the ToF of the paths. When this function if called, the vector is appended by one more entry which corresponds to the ToF of the new path
% angleFromMusicPrev: (number of paths calculated so far)-dimensional vector containing the AoA of the paths. When this function if called, the vector is appended by one more entry which corresponds to the AoA of the new path
% deltaFromMusicPrev: (number of paths calculated so far)-dimensional vector containing the AoD of the paths. When this function if called, the vector is appended by one more entry which corresponds to the AoD of the new path
% elevationFromMusicPrev: same as above but for AoD elevation angle in degrees
% SubCarrInd: indices fo WiFi subcarriers at which CSI is calculated
% deltaGridValue: MUSIC spectrum value is calculated for different possible values of the propagation path parameters. This vector contains the values of AoD parameters for which MUSIC algorithm is calculated.
% u_sGridValue: same as above for AoA values
% delayGridValue: same as above for ToF values
% elevationGridValue: same as above for AoD elevation degrees
% T: number of transmit antennas
% c: speed of light
% do_second_iter: use 0
% seconditerGridPts: not valid
% useNoise: use 0
% dTx: distance of transmitter antennas
% % % OUTPUT:
% delayFromMusic: (number of paths calculated so far by end of this function execution)-dimensional vector containing the AoA of the paths. 
% angleFromMusic: same as above for AoA
% deltaFromMusic: same as above for AoD
% elevationFromMusic: same as above for elevation angle of AoD
% maxCorr: maximum value of MUSIC spectrum                                                
% music_spectrum_plot: MUSIC spectrum                                                
                                                
                                                
                                                
                                                
maxCorr = [];
if ~useNoise
    if ~isempty(delayFromMusicPrev)
        Ahat = [];
        for compNo = 1:length(delayFromMusicPrev)
            u_s = (d*fc/c)*sind(angleFromMusicPrev(compNo));
            Ahat(:,compNo) = circularGridSampleBackscatter3D(fc, 0, 0, K, u_s, c, SubCarrInd(1:L), fgap, delayFromMusicPrev(compNo), dTx, deltaFromMusicPrev(compNo), elevationFromMusicPrev(compNo) );
        end
        PerpAhat = eye(size(Ahat,1)) - Ahat*pinv(Ahat'*Ahat)*Ahat';
    else
        PerpAhat = eye(size(Qs,1),'like',2+1i);
    end
    Qs = PerpAhat*Qs;
    
    chirp_z = 0;
    doHadamard = 1;
    if chirp_z
        QsConj = conj(Qs);
        QsReshape = reshape(QsConj,L,K,[],size(QsConj,2));
        
        deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);
        aoaSteeringInvMat = exp(1i*2*pi*dTx*(fc/c)*[zeros(size(deltaConsider)); -cosd(deltaConsider + 60); -cosd(deltaConsider)]);
        u_sConsider = (d*fc/c)*sind( GridStart(2) + (0:GridPts(2)-1)*GridSpacing(2) );
        aoaSteeringMat = exp(-1i*2*pi*((-(K-1)/2:(K-1)/2)')*u_sConsider);
        music_spec_num = musicChirpSpectrum(QsReshape, aoaSteeringInvMat, aoaSteeringMat, GridPts, GridSpacing, GridStart,...
                                    SubCarrInd, fgap, d, fc, c);
        
        PerpAhatReshape = reshape(PerpAhat.',L,K,[],size(PerpAhat.',2));
        music_spec_den = musicChirpSpectrum(PerpAhatReshape, aoaSteeringInvMat, aoaSteeringMat, GridPts, GridSpacing, GridStart,...
                                    SubCarrInd, fgap, d, fc, c);
     
        music_spectrum = abs(music_spec_num./music_spec_den);
    end
    if ~chirp_z
        if ~doHadamard
            Ps = (Qs*Qs');

            % % achieving the music spectrum without for loop by only matrix
            % % computations
            % delayConsider = GridStart(1) + (0:GridPts(1)-1)*GridSpacing(1);
            % u_sConsider = (d*fc/c)*sind( GridStart(2) + (0:GridPts(2)-1)*GridSpacing(2) );
            % deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);
            % 
            % delaySteeringMat = exp(-1i*2*pi*(SubCarrInd(1:L).')*fgap*delayConsider);
            % aoaSteeringMat = exp(-1i*2*pi*((-(K-1)/2:(K-1)/2)')*u_sConsider);
            % aoaSteeringInvMat = exp(1i*2*pi*dTx*(fc/c)*[zeros(size(deltaConsider)); -cosd(deltaConsider + 60); -cosd(deltaConsider)]);
            % 
            % thetaTauMat = kron(aoaSteeringMat, delaySteeringMat);
            % thetaTauDeltaMat = kron(aoaSteeringInvMat, thetaTauMat);
            % 
            % % we have now created steering vectors where each column is a steering
            % % vector.
            % PerpAhat_a = PerpAhat*thetaTauDeltaMat; % modifying each steering vector. These are the a1's.
            % PsA = Ps*PerpAhat_a;
            % music_spec_num = sum((PerpAhat_a').*(PsA.'),2);            % now conjugate of each a1 should be dot producted with columns of PsA
            % music_spec_den = sum((PerpAhat_a').*(PerpAhat_a.'),2); % now conj of each a1 should be dot producted with a1
            % music_spectrum = abs(music_spec_num./music_spec_den);

            % previous way of computing MUSIC spectrum by for looping over all
            % possible grid values
            numGridPts = prod(GridPts);
            music_spectrum = zeros(numGridPts,1);
            % Use PARFOR
            if ~isempty(aTot)
                aTotFull = 1;
            else
                aTotFull = 0;
            end
            for loopVar = 1:numGridPts
                if aTotFull
                    aTotTmp = aTot(:,loopVar); % use this if you have precomputed the grid
                else
                    aTotTmp = circularGridSampleBackscatter3D(fc, T, deltaGridValue(loopVar), K, u_sGridValue(loopVar), c, SubCarrInd(1:L), fgap, delayGridValue(loopVar), dTx, deltaGridValue(loopVar), elevationGridValue(loopVar) );
                end
                a1 = PerpAhat*aTotTmp;
                music_spectrum(loopVar) = (a1'*Ps*a1)/(a1'*a1); 
            end
            music_spectrum = abs(music_spectrum);
        end
        if doHadamard
            A = PerpAhat * aTot;
            QA = Qs'*A;
            music_spectrum = norms(QA,2,1)./norms(A,2,1);
            music_spectrum = music_spectrum.^2;
        end
    end
    
else
    %     % delay forms first dimension, angle second and displacement third
    %     delayConsider = GridStart(1) + (0:GridPts(1)-1)*GridSpacing(1);
    %     u_sConsider = (d*fc/c)*sind( GridStart(2) + (0:GridPts(2)-1)*GridSpacing(2) );
    %     deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);
    %     % delayConsider = unique(delayGridValue, 'stable');
    %     % u_sConsider = unique(u_sGridValue, 'stable');
    %     % deltaConsider = unique(deltaGridValue, 'stable');
    %     
    %     delaySteeringMat = exp(-1i*2*pi*(SubCarrInd(1:L).')*fgap*delayConsider);
    %     aoaSteeringMat = exp(-1i*2*pi*((-(K-1)/2:(K-1)/2)')*u_sConsider);
    %     aoaSteeringInvMat = exp(1i*2*pi*dTx*(fc/c)*[zeros(size(deltaConsider)); -cosd(deltaConsider + 60); -cosd(deltaConsider)]);
    %     
    %     thetaTauMat = kron(aoaSteeringMat, delaySteeringMat);
    %     thetaTauDeltaMat = kron(aoaSteeringInvMat, thetaTauMat);
    %     
    %     PnA = (Qn*Qn')*thetaTauDeltaMat;
    %     thetaTauDeltaMatTrans = thetaTauDeltaMat';
    %     % we want only the diagonal terms of thetaTauDeltaMatTrans*PnA
    %     music_spectrum = sum(thetaTauDeltaMatTrans.*(PnA.'),2);
    % 
    %     music_spectrum = 1./abs(music_spectrum);
    %     % figure; imagesc(reshape(squeeze(music_spectrum), GridPts([1 3])))
    %     % figure; imagesc(reshape((music_spectrum), GridPts([1 2])))
    
    A = aTot;
    QA = Qn'*A;
    music_spectrum = norms(A,2,1)./norms(QA,2,1);
    music_spectrum = music_spectrum.^2;
end

%% Finding maximum of Music specturm and the corresponding Grid values
if ~useNoise
    [maxCorr,loopVarMax] = max(music_spectrum); % returns the first occurence of the maximum
    % if doGPU
    %     maxCorr = gather(maxCorr);
    %     loopVarMax = gather(loopVarMax);
    % end
    [delay_idx,angle_idx, delta_idx, elevation_idx] = ind2sub(GridPts,loopVarMax);
    delayFromMusic = GridStart(1) + (delay_idx-1)*GridSpacing(1);
    angleFromMusic = GridStart(2) + (angle_idx-1)*GridSpacing(2);
    deltaFromMusic = GridStart(3) + (delta_idx-1)*GridSpacing(3);
    elevationFromMusic = GridStart(4) + (elevation_idx-1)*GridSpacing(4);
else
    % assuimg that ToF does not come into picture, only GridPts[2 3 4] matter
    maxPeaks = Inf;
    music_spectrum = reshape(music_spectrum,GridPts(2),GridPts(3),GridPts(4));
    BW = imregionalmax(music_spectrum);
    [angle_idx, delta_idx, elevation_idx] = ind2sub(size(BW), find(BW));
    topPeakIndices = topk(music_spectrum(BW), min(maxPeaks, length(find(BW(:)))));
    delayFromMusic = GridStart(1)*ones(length(topPeakIndices),1);
    angleFromMusic = GridStart(2) + (angle_idx(topPeakIndices)-1)*GridSpacing(2);
    deltaFromMusic = GridStart(3) + (delta_idx(topPeakIndices)-1)*GridSpacing(3);
    elevationFromMusic = GridStart(4) + (elevation_idx(topPeakIndices)-1)*GridSpacing(4);
    musicLocalMaxInd = find(BW);
    maxCorr = music_spectrum(musicLocalMaxInd(topPeakIndices));
end

music_spectrum = reshape(music_spectrum,GridPts(1),GridPts(2),GridPts(3),GridPts(4));
music_spectrum_plot = music_spectrum;

%% for second iteration of MUSIC

if do_second_iter  
    
    delayFirstIter = delayFromMusic;
    angleFirstIter = angleFromMusic;
    deltaFirstIter = deltaFromMusic;
    elevationFirstIter = elevationFromMusic;
    
    for iComp = 1:length(delayFirstIter)
        GridPts = seconditerGridPts;
        delayRange = [delayFirstIter(iComp)-GridSpacing(1) delayFirstIter(iComp)+GridSpacing(1)];
        angleRange = [angleFirstIter(iComp)-GridSpacing(2) angleFirstIter(iComp)+GridSpacing(2)];
        deltaRange = [deltaFirstIter(iComp)-GridSpacing(3) deltaFirstIter(iComp)+GridSpacing(3)];
        elevationRange = [elevationFirstIter(iComp)-GridSpacing(4) elevationFirstIter(iComp)+GridSpacing(4)];
        [GridStart, GridSpacing, delayGridValue, u_sGridValue, deltaGridValue, elevationGridValue] = gridParamsBackscatter3D(GridPts, angleRange, deltaRange, d, fc, c, delayRange, elevationRange);

        if ~useNoise
            %numGridPts = prod(GridPts);
            if chirp_z
                QsConj = conj(Qs);
                QsReshape = reshape(QsConj,L,K,[],size(QsConj,2));

                deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);
                aoaSteeringInvMat = exp(1i*2*pi*dTx*(fc/c)*[zeros(size(deltaConsider)); -cosd(deltaConsider + 60); -cosd(deltaConsider)]);
                u_sConsider = (d*fc/c)*sind( GridStart(2) + (0:GridPts(2)-1)*GridSpacing(2) );
                aoaSteeringMat = exp(-1i*2*pi*((-(K-1)/2:(K-1)/2)')*u_sConsider);
                music_spec_num = musicChirpSpectrum(QsReshape, aoaSteeringInvMat, aoaSteeringMat, GridPts, GridSpacing, GridStart,...
                                            SubCarrInd, fgap, d, fc, c);

                PerpAhatReshape = reshape(PerpAhat.',L,K,[],size(PerpAhat.',2));
                music_spec_den = musicChirpSpectrum(PerpAhatReshape, aoaSteeringInvMat, aoaSteeringMat, GridPts, GridSpacing, GridStart,...
                                            SubCarrInd, fgap, d, fc, c);

                music_spectrum = abs(music_spec_num./music_spec_den);
            end
            if ~chirp_z
                % % achieving the music spectrum without for loop by only matrix
                % % computations
                % delayConsider = GridStart(1) + (0:GridPts(1)-1)*GridSpacing(1);
                % u_sConsider = (d*fc/c)*sind( GridStart(2) + (0:GridPts(2)-1)*GridSpacing(2) );
                % deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);
                % 
                % delaySteeringMat = exp(-1i*2*pi*(SubCarrInd(1:L).')*fgap*delayConsider);
                % aoaSteeringMat = exp(-1i*2*pi*((-(K-1)/2:(K-1)/2)')*u_sConsider);
                % aoaSteeringInvMat = exp(1i*2*pi*dTx*(fc/c)*[zeros(size(deltaConsider)); -cosd(deltaConsider + 60); -cosd(deltaConsider)]);
                % 
                % thetaTauMat = kron(aoaSteeringMat, delaySteeringMat);
                % thetaTauDeltaMat = kron(aoaSteeringInvMat, thetaTauMat);
                % 
                % % we have now created steering vectors where each column is a steering
                % % vector.
                % PerpAhat_a = PerpAhat*thetaTauDeltaMat; % modifying each steering vector. These are the a1's.
                % PsA = Ps*PerpAhat_a;
                % music_spec_num = sum((PerpAhat_a').*(PsA.'),2);            % now conjugate of each a1 should be dot producted with columns of PsA
                % music_spec_den = sum((PerpAhat_a').*(PerpAhat_a.'),2); % now conj of each a1 should be dot producted with a1
                % music_spectrum = abs(music_spec_num./music_spec_den);

                % %         Previous way of doing it through far loop
                numGridPts = prod(GridPts);
                music_spectrum = zeros(numGridPts,1);
                for loopVar = 1:numGridPts
                    aTotTmp = circularGridSampleBackscatter3D(fc, T, deltaGridValue(loopVar), K, u_sGridValue(loopVar), c, SubCarrInd(1:L), fgap, delayGridValue(loopVar), dTx, deltaGridValue(loopVar), elevationGridValue(loopVar) );
                    a1 = PerpAhat*aTotTmp;
                    music_spectrum(loopVar) = (a1'*Ps*a1)/(a1'*a1); 
                end
                music_spectrum = abs(music_spectrum);
            end
        else
            % delay forms first dimension, angle second and displacement third
            delayConsider = GridStart(1) + (0:GridPts(1)-1)*GridSpacing(1);
            u_sConsider = (d*fc/c)*sind( GridStart(2) + (0:GridPts(2)-1)*GridSpacing(2) );
            deltaConsider = GridStart(3) + (0:GridPts(3)-1)*GridSpacing(3);
            % delayConsider = unique(delayGridValue, 'stable');
            % u_sConsider = unique(u_sGridValue, 'stable');
            % deltaConsider = unique(deltaGridValue, 'stable');

            delaySteeringMat = exp(-1i*2*pi*(SubCarrInd(1:L).')*fgap*delayConsider);
            aoaSteeringMat = exp(-1i*2*pi*((-(K-1)/2:(K-1)/2)')*u_sConsider);
            aoaSteeringInvMat = exp(1i*2*pi*dTx*(fc/c)*[zeros(size(deltaConsider)); -cosd(deltaConsider + 60); -cosd(deltaConsider)]);

            thetaTauMat = kron(aoaSteeringMat, delaySteeringMat);
            thetaTauDeltaMat = kron(aoaSteeringInvMat, thetaTauMat);

            PnA = (Qn*Qn')*thetaTauDeltaMat;
            thetaTauDeltaMatTrans = thetaTauDeltaMat';
            % we want only the diagonal terms of thetaTauDeltaMatTrans*PnA
            music_spectrum = sum(thetaTauDeltaMatTrans.*(PnA.'),2);

            music_spectrum = 1./abs(music_spectrum); 
        end

        if ~useNoise
            [maxCorrTmp,loopVarMaxTmp] = max(music_spectrum); 
            maxCorr = maxCorrTmp(1);
            loopVarMax = loopVarMaxTmp(1);
            [delay_idx,angle_idx, delta_idx, elevation_idx] = ind2sub(GridPts,loopVarMax); 
            delayFromMusic = GridStart(1) + (delay_idx-1)*GridSpacing(1);
            angleFromMusic = GridStart(2) + (angle_idx-1)*GridSpacing(2);
            deltaFromMusic = GridStart(3) + (delta_idx-1)*GridSpacing(3);
            elevationFromMusic = GridStart(4) + (elevation_idx-1)*GridSpacing(4);
        else
            maxCorr = max(music_spectrum);
            music_spectrum = reshape(music_spectrum,GridPts(1),GridPts(2),GridPts(3));
            BW = imregionalmax(music_spectrum);
            [delay_idx, angle_idx, delta_idx] = ind2sub(size(BW), find(BW));
            topPeakIndices = topk(music_spectrum(BW), 1);
            delayFromMusic(iComp) = GridStart(1) + (delay_idx(topPeakIndices)-1)*GridSpacing(1);
            angleFromMusic(iComp) = GridStart(2) + (angle_idx(topPeakIndices)-1)*GridSpacing(2);
            deltaFromMusic(iComp) = GridStart(3) + (delta_idx(topPeakIndices)-1)*GridSpacing(3);
        end
    end
end

