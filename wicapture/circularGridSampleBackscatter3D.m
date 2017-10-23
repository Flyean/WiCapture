%% produces steering vector for a given set of propagation path parameters
function steeringVec = circularGridSampleBackscatter3D(fc, T, deltak, M, u_s, c, SubCarrInd, fgap, delay, dTx, dodRay, dodRayElevation )
% % % INPUTS:
% fc: center frequency
% T: number of time instants with equal displacements
% deltak: dispacement along the thetak direction of the transmitter
% M: number of receive antennas
% u_s: (dstance between adjacent receive antennas)*sin(AoA of different paths)/(wavelength)
% c: speed of light
% SubCarrInd: WiFi subcarrier indices at which CSI is measured
% fgap: frequency gap between adjacent WiFi subcarriers
% delay: ToF of all the paths
% dTx: distance between adjacent antennas at the transmitter
% dodRay: AoD-azimuth(angle of departure) of paths in degrees
% dodRayElevation: AoD-elevation(angle of departure) of paths in degrees
% % % OUTPUT:
% steeringVec: steering vector

aoaSteering = exp(-1i*2*pi*u_s*(0:M-1)');
% aoaSteering = exp(-1i*2*pi*u_s*(-(M-1)/2:(M-1)/2)'); % finding AoA steering vector using symmetric indices -old times
dodSteering = exp(1i*2*pi*dTx*(fc/c)* sind(dodRayElevation) * [0 -cosd((dodRay + 60)) -cosd((dodRay))]');
%aoaSteeringInv = exp(-1i*2*pi*fc*(0:(T-1))'*deltak/c);
delaySteering = exp(-1i*2*pi*(SubCarrInd(:))*fgap*delay);

steeringVecDelayAoATmp = delaySteering*aoaSteering.';
steeringVecDelayAoA = steeringVecDelayAoATmp(:);
steeringVecTmp = steeringVecDelayAoA*dodSteering.';
steeringVec = steeringVecTmp(:);
