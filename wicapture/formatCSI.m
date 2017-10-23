%% formatting the input CSI matrix
function X = formatCSI(csi_trace, N, M, Ttot, K, L, T)
% % % INPUT: 
% N: number of such subcarriers
% M: number of receive antennas
% Ttot: number of transmit antennas
% K: number of receive antennas in 'smoothed CSI'
% L: number of subcarriers in 'smoothed CSI'
% T: number of transmit antennas in smoothed CSI
% csi_trace : (N*M*Ttot)x1 vector containing the csi trace
% % % OUTPUT:
% X is smoothed CSI matrix that we input to MUSIC 

%% reshaping CSI to form smoothed matrix
% X is the matrix which we input to MUSIC
X = [];
X3D = [];
csiTracePerPkt = reshape(csi_trace, N, M, Ttot);
for iT = 1:Ttot
    ysenarr = csiTracePerPkt(:,:,iT);
    [N,M] = size(ysenarr);
    D = [];
    for m=1:M %mth antenna
        D{m} = hankel(ysenarr(1:L,m),ysenarr(L:N,m)); 
    end
    De = zeros(K*L, (M - K + 1)*(N - L + 1));
    for start_idx = 1:K
        tmp =  cat(2,D{start_idx:start_idx+M-K});
        De((start_idx-1)*L+(1:L),:) = tmp;
    end
    X = [X; De];
    X3D{iT} = De;
end
% % % if 3D MUSIC needs to be performed
X = zeros(K*L*T, (M - K + 1)*(N - L + 1)*(Ttot - T + 1));
subarraySizeAntSubcarr = L*K;
for start_idx = 1:T
    tmp =  cat(2,X3D{start_idx:start_idx+Ttot-T});
    X((start_idx-1)*subarraySizeAntSubcarr + (1:subarraySizeAntSubcarr),:) = tmp;
end
