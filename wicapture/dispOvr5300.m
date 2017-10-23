function [trajEst, dispEst_fromalpha, A_dod_comb_allRxSave, compMinVar_dispDiffMusic_comb_allRxSave] = dispOvr5300(aoaDodOnly,expStruct, nSegments, nRx, txOrientation)
% % % INPUT:
% aoaDodOnly: (number of receivers x number of WiFi packets) dimensional cell which contains the propgation path parameters
% expStruct: a structure containing experiment parameters like center frequency. Look at wicaptureAlgo.m for more information of this structure.
% nSegments: number of WiFi packets. Position is estimated at each of the packets.
% nRx: number of receivers
% txOrientation: orientation of the transmitter. just feed 0 for this parameter.
% % % OUTPUT:
% trajEst: estimate of the WICapture trajectory
% dispEst_fromalpha: estimate of dispacement between consecutive packets
% A_dod_comb_allRxSave: matrix in equation (12) in the paper formed from the AoD of different paths
% compMinVar_dispDiffMusic_comb_allRxSave: change in the phase of complex attenuation of the paths 

v2struct(expStruct);

trajTrue = zeros(nSegments,2);
% displacements from trajectory
dispTraj = (diff(trajTrue)).';

A_dod_comb_allRx = [];
compMinVar_dispDiffMusic_comb_allRx =[];
weightVec_comb_allRx = [];
compMinVar_id = cell(nRx,1);
topComps = cell(nRx,1);

for iRx = 1:nRx
    txGap = 1;
    %% differences of phases across locations from complex amplitudes at each subcarrier
    paramEachSubcarr = cat(1, aoaDodOnly{iRx,:});
    ampEachSubcarr = paramEachSubcarr(:,alphaInd);
    % create (NxnSegmentsxL) matrix of complex amplitudes - ampEachComp
    Lest = size(aoaDodOnly{iRx,1},1); % estimate of the number of components
    ampEachCompTmp = reshape(ampEachSubcarr, Lest, N, nSegments);
    ampEachComp = permute(ampEachCompTmp, [2 3 1]);
    phase_allSubcarr = angle(exp(1i*(angle(ampEachComp))));
    phaseDiffLoc_allSubcarr = angle(exp(1i*(diff(angle(ampEachComp),1,2)))); % for all subcarriers, find the difference in the phases at successive locations

    phaseLoc = zeros(Lest, size(phaseDiffLoc_allSubcarr,2)); % average phase shift for successive locations, averaged over all subcarriers
    phaseDiffLoc = zeros(Lest, nSegments-txGap);
    phaseDiffLocVar = phaseDiffLoc;
    meanAmpLoc = phaseDiffLoc;
    ampDiffLocVar = phaseDiffLoc;
    phaseLinMeas = phaseDiffLoc;
    % % % % Commenting to improve speed
    % for iLoc = 1:nSegments
    %     for iComp = 1:Lest
    %         pCoeff = polyfit(SubCarrInd(:),unwrap(phase_allSubcarr(:,iLoc,iComp)),1);
    %         phaseLoc(iComp, iLoc) = pCoeff(end); % finding avg phase change across locations
    %         phaseLinMeas(iComp, iLoc) = norm(unwrap(phase_allSubcarr(:,iLoc,iComp)) - pCoeff(1)*SubCarrInd(:) - pCoeff(2))/sqrt(N);
    %     end
    % end
    
    for iLoc = 1:nSegments-txGap
        for iComp = 1:Lest  
            % % % % Commenting to improve speed
            % stoCfoAngle = unwrap(angle(exp(1i*(phase_allSubcarr(:,iLoc+txGap,iComp) - phase_allSubcarr(:,iLoc,iComp)))));
            % pCoeffDiff = polyfit(SubCarrInd(:),stoCfoAngle,1);
            % phaseDiffLoc(iComp, iLoc) = angle(exp(1i*(pCoeffDiff(end)))); % unwrapping changes the phase in multiples of 2*pi
            % phaseDiffLocVar(iComp, iLoc) = norm(stoCfoAngle - pCoeffDiff(1)*SubCarrInd(:) - pCoeffDiff(2))/sqrt(N);
            
            phaseDiffLoc(iComp, iLoc) = mean(angle(exp(1i*(phase_allSubcarr(:,iLoc+txGap,iComp) - phase_allSubcarr(:,iLoc,iComp)))));
            % phaseDiffLoc(iComp, iLoc)
            % phaseDiffLoc(iComp, iLoc) = mean(phase_allSubcarr(:,iLoc+txGap,iComp) - phase_allSubcarr(:,iLoc,iComp));
            % phaseDiffLoc(iComp, iLoc)
            % phaseDiffLocVar(iComp, iLoc) = std(phase_allSubcarr(:,iLoc+txGap,iComp) - phase_allSubcarr(:,iLoc,iComp));
            % figure(2); 
            % plot(phase_allSubcarr(:,iLoc+txGap,iComp) - phase_allSubcarr(:,iLoc,iComp)); hold on;
            % plot(phaseDiffLoc_allSubcarr(:,iLoc,iComp))
            % plot(stoCfoAngle); hold off;
                        
            meanAmpLoc(iComp, iLoc) = norm(ampEachComp(:,iLoc,iComp));%abs(mean(ampEachComp(:,iLoc+txGap,iComp)./ampEachComp(:,iLoc,iComp)));
            ampDiffLocVar(iComp, iLoc) = std(ampEachComp(:,iLoc+txGap,iComp)./ampEachComp(:,iLoc,iComp));

            % figure(11);
            % plot(phase_allSubcarr(:,iLoc+txGap,iComp))
            % hold on
            % plot(phase_allSubcarr(:,iLoc,iComp))
            % plot(angle(exp(1i*(stoCfoAngle))))
            % % plot(-(- pCoeffDiff(1)*SubCarrInd(:) - pCoeffDiff(2)))
            % % plot(phaseDiffLoc_allSubcarr(:,iLoc,iComp));
            % hold off;
            % ylim(pi*[-1 1])
            % title(sprintf('comp %d loc %d', iComp, iLoc))
            
        end
    end
    
    %     % debug
    %     for iLoc = 1:nSegments-txGap
    %         for iComp = 1:Lest    
    % %             phaseLoc(iComp, iLoc+txGap) - phaseLoc(iComp, iLoc)
    % %             phaseDiffLoc(iComp, iLoc)
    %             figure(100);
    %             plot(phase_allSubcarr(:,iLoc,iComp)); hold on;
    %             plot(phase_allSubcarr(:,iLoc+txGap,iComp)); hold on;
    %             plot(phaseDiffLoc_allSubcarr(:,iLoc,iComp)); hold off;
    % %             sprintf('iComp %d iLoc %d', iComp, iLoc)
    %         end
    %     end
    
    % finding difference in phase for all multipath
    compMinVar_dispDiffMusic = [];
    for iExp = 1:nSegments-txGap
        compMinVar_dispDiffMusic(:, iExp) = -phaseDiffLoc(:, iExp)/(2*pi*fc/c);   
    end

    % if you have more than two paths, use the two paths with largest amplitude
    topComps{iRx} = min(2,Lest);
    compMinVar_id{iRx} = topk( mean(meanAmpLoc,2)./std(meanAmpLoc,[],2),topComps{iRx});
    %compMinVar_id{iRx} = topk(mean(meanAmpLoc,2),topComps{iRx});
    % compMinVar_id = topk(1./sum(phaseDiffLocVar,2),topComps);
    topComps{iRx} = length(compMinVar_id{iRx});

    % taking only the topComps number of multipaths for displacement calculation
    compMinVar_dispDiffMusic = compMinVar_dispDiffMusic(compMinVar_id{iRx},:);
    weightVec = (mean(meanAmpLoc(compMinVar_id{iRx},:),2));

    % subtracting the rows of difference of phase measurements corresponding to different multipath components
    compMinVar_dispDiffMusic_comb = [];
    weightVec_comb = [];
    for iComp = 1:topComps{iRx}
        for jComp = iComp:topComps{iRx}
            if iComp~=jComp
                compMinVar_dispDiffMusic_comb(end+1,:) = angle(exp(1i*(2*pi*fc/c)*(compMinVar_dispDiffMusic(iComp,:) - compMinVar_dispDiffMusic(jComp,:))))/(2*pi*fc/c);
                weightVec_comb(end+1,:) = weightVec(iComp)*weightVec(jComp);                
            end
        end
    end
    compMinVar_dispDiffMusic_comb_allRx = [compMinVar_dispDiffMusic_comb_allRx; compMinVar_dispDiffMusic_comb];
    weightVec_comb_allRx = [weightVec_comb_allRx; weightVec_comb];

end


% using difference of measurements with weights
dispEst_fromalpha = [];
A_dod_comb_allRxSave = cell(nSegments-txGap,1);
for iExp = 1:nSegments-txGap
    A_dod_comb_allRx = [];
    for iRx = 1:nRx        
        dod_rayEst = aoaDodOnly{iRx,iExp}(compMinVar_id{iRx},dodInd) + txOrientation;
        A_dod = -[cosd(dod_rayEst) sind(dod_rayEst)];
        
        % subtracting the rows of A_dod corresponding to different multipath components
        A_dod_comb = [];
        for iComp = 1:topComps{iRx}
            for jComp = iComp:topComps{iRx}
                if iComp~=jComp
                    A_dod_comb(end+1,:) = A_dod(iComp,:) - A_dod(jComp,:); 
                end
            end
        end
        A_dod_comb_allRx = [A_dod_comb_allRx; A_dod_comb];
    end
    %     % % % % using lower and upper bounds on data
    %     options = optimset('Display','none','MaxIter',500);
    %     lb = -[5 5]*1e-3;
    %     dispEst_fromalpha(:,iExp) = lsqlin(A_dod_comb_allRx,compMinVar_dispDiffMusic_comb_allRx(:, iExp),[],[],[],[],lb,-lb,[],options);
    % % % alternative for CVX
    dispEst_fromalpha(:,iExp) = lscov(A_dod_comb_allRx, compMinVar_dispDiffMusic_comb_allRx(:, iExp), weightVec_comb_allRx);
    % % % % CVX method
    % cvx_begin quiet
    %     variable dispTmp(2)
    %     minimize(norm(weightVec_comb_allRx.*(A_dod_comb_allRx*dispTmp - compMinVar_dispDiffMusic_comb_allRx(:, iExp))))
    %     subject to
    %         norm(dispTmp) <= c/fc/4
    % cvx_end
    % dispEst_fromalpha(:,iExp) = dispTmp;
    % % % % using difference of linear equations from multiple APs
    % dispEst_fromalpha(:,iExp) = lscov(A_dod_comb_allRx, compMinVar_dispDiffMusic_comb_allRx(:, iExp), weightVec_comb_allRx);
    % % % % naive solving of differences of linear equations
    % dispEst_fromalpha(:,iExp) = lscov(A_dod_comb, compMinVar_dispDiffMusic_comb(:, iExp));
    % % % % solving linear equations with weights 
    % dispEst_fromalpha(:,iExp) = lscov(A_dod, compMinVar_dispDiffMusic(:, iExp), (mean(meanAmpLoc(compMinVar_id,:),2)));
    % dispEst_fromalpha(:,iExp) = linsolve(A_dod, compMinVar_dispDiffMusic(:, iExp));
    
    A_dod_comb_allRxSave{iExp} = A_dod_comb_allRx;
    compMinVar_dispDiffMusic_comb_allRxSave = compMinVar_dispDiffMusic_comb_allRx;
end

% % trajectory estimation
trajEst = 1e3*cumsum([zeros(2,1) dispEst_fromalpha].');
    
    