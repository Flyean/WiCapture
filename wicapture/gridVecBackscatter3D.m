%% creates set of steering vectors. MUSIC spectrum is calculated at each of these steering vectors.
function [aTot,GridStart,GridSpacing, delayGridValue, u_sGridValue, deltaGridValue, elevationGridValue] = gridVecBackscatter3D(deltaRange, M, T, d, fc, c, do_second_iter, delayRange, SubCarrInd, N, fgap, GridPts, MaxAngle, MinAngle, elevationRange, generateAtot)
% % % INPUT:
% deltaRange: 1x2 dimensional vector with the minimum and maximum values of the AoD azimuth angle in degrees for which MUSIC spectrum needs to be calculated.
% M: number of receive antennas
% T: number of transmit antennas
% fc: center frequency
% c: speed of light
% do_second_iter: use 0
% delayRange: 1x2 dimensional vector with the minimum and maximum values of the ToF in s for which MUSIC spectrum needs to be calculated.
% SubCarrInd: indices of subcarriers where CSI is provided
% N: number of such subcarriers
% fgap: frequency gap between successive subcarriers
% GridPts: 1x4 vector: number of points along each dimensions of [ToF AoA AoD AoD_elevation] parameters at which MUSIC spectrum is estimated.
% MaxAngle: maximum value of AoA in degrees for which MUSIC spectrum needs to be calculated.
% MinAngle: minimum value of AoA in degrees for which MUSIC spectrum needs to be calculated.
% elevationRange: 1x2 dimensional vector with the minimum and maximum values of the AoD elevation angle in degrees for which MUSIC spectrum needs to be calculated.
% generateAtot: 0 if you want to use matrix of steering vectors from file, 1 if you want to recalculate the set of steering vectors, 2 to recalculate and save to a file.
% % % OUTPUT:
% aTot: matrix with each column being a steering vector
% GridStart: 1x4 vector: smallest value of each of [ToF AoA AoD AoD_elevation] parameters at which MUSIC spectrum is estimated.
% GridSpacing: 1x4 vector: difference between successive grid values of [ToF AoA AoD AoD_elevation] parameters at which MUSIC spectrum is estimated.
% deltaGridValue: MUSIC spectrum value is calculated for different possible values of the propagation path parameters. This vector contains the values of AoD parameters for which MUSIC algorithm is calculated.
% u_sGridValue: same as above for AoA values
% delayGridValue: same as above for ToF values
% elevationGridValue: same as above for AoD elevation degrees


% do_second_iter = 0;

% Define angles at which MUSIC ?spectrum? will be computed
GridStart = [delayRange(1) MinAngle deltaRange(1) elevationRange(1)];
GridStop = [delayRange(2) MaxAngle deltaRange(2) elevationRange(2)];

% GridSpacing = (GridStop-GridStart)./(GridPts-ones(size(GridPts)));
GridSpacing = (GridStop-GridStart)./max(1, GridPts-ones(size(GridPts)));

angleStart = GridStart(2);
delayStart = GridStart(1);
deltaStart = GridStart(3);
elevationStart = GridStart(4);
angleDiff = GridSpacing(2);
delayDiff = GridSpacing(1);
deltaDiff = GridSpacing(3);
elevationDiff = GridSpacing(4);
% siz = GridPts;
numGridPoints = prod(GridPts);

[delayIndices,angleIndices,deltaIndices,elevationIndices] = ind2sub(GridPts,1:numGridPoints);
delayGridValue = delayStart + (delayIndices-1)*delayDiff;
angleGridValue = (angleStart + (angleIndices-1)*angleDiff)*pi/180;
u_sGridValue = (d*fc/c)*sin(angleGridValue);
deltaGridValue = deltaStart + (deltaIndices-1)*deltaDiff;
elevationGridValue = elevationStart + (elevationIndices-1)*elevationDiff;

aTot = [];
% % % If you want the aTot to be nonempty, uncomment below
%generateAtot = 0;
if ~generateAtot
    load('aTotSave');
else
    aTot = zeros(M*T,numGridPoints);
    for loopVar = 1:numGridPoints
        aTot(:,loopVar) = circularGridSampleBackscatter3D(fc, T, deltaGridValue(loopVar), M, u_sGridValue(loopVar),...
            c, 0, fgap, delayGridValue(loopVar), d, deltaGridValue(loopVar), elevationGridValue(loopVar) );
    end
    if generateAtot == 2
        save('aTotSave','aTot')
    end
end
