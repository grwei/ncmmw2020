%% Basic information
% The verification of the following paper(Section III):
% M. K. Sampath, "On addressing the practical issues in the extraction of RLGC parameters for lossy multiconductor transmission lines using S-parameter models," 2008 IEEE-EPEP Electrical Performance of Electronic Packaging, San Jose, CA, 2008, pp. 259-262.
% Author Name: Wei Guorui
% Instructor Name: Xia Bin
% Created at: 2019-10-24 09:52
% Last modified: 2019-10-24 
clc; clear; close all;
addpath(genpath('../code'))

%% Port numbering
% 	           PORT NUMBERING (Classic style)
% 	    PORT#                               PORT#
% 	    Z1 = Z3 = 50.00                     Z2 = Z4 = 50.00
% 	           --------------------------
% 	          |                          |
% 	      1 --|                          |-- N+1
% 	          |                          |
% 	     ...--|                          |-- ...
% 	          |                          |
% 	    N-1 --|                          |-- 2N-1
% 	          |                          |
% 	      N --|                          |-- 2N
% 	          |                          |
% 	           --------------------------

%% Initialize
LL_1 = 0.0254; % line length(meters).
LL_2 = 0.0224; 
LL = LL_1 - LL_2;
filename_TL1 = 'data/SI/d4line_25_4_mm_201910242117.s4p';
filename_TL2 = 'data/SI/d4line_22_4_mm_201910242117.s4p';
SE8P_Data_1 = read(rfdata.data,filename_TL1); % SingleEnded8PortData
SE8P_Data_2 = read(rfdata.data,filename_TL2); % SingleEnded8PortData
Z0 = SE8P_Data_1.Z0;
ZL = SE8P_Data_1.ZL;
ZS = SE8P_Data_1.ZS;
SE8P_Data_1.S_Parameters = snp2smp(SE8P_Data_1.S_Parameters,...
    Z0,[1 2 3 4]);
SE8P_Data_2.S_Parameters = snp2smp(SE8P_Data_2.S_Parameters,...
    Z0,[1 2 3 4]);
freq = SE8P_Data_1.Freq;
freqpts = length(freq);
numLines = size(SE8P_Data_1.S_Parameters,1)/2;
s_params_1 = sparameters(SE8P_Data_1);
s_params_2 = sparameters(SE8P_Data_2);

%% Convert 2N¡Á2N S matrix to 2N¡Á2N T Matrix(ABCD) and Z matrix
T_mat_1 = s2t_mod(SE8P_Data_1.S_Parameters, Z0); % notice the difference in "s2rlgc": S -> T.
T_mat_2 = s2t_mod(SE8P_Data_2.S_Parameters, Z0); % notice the difference in "s2rlgc": S -> T.

%% Calculate the M matrix
X_mat = T_mat_1 + T_mat_2;
Y_mat = T_mat_1 - T_mat_2;
M_mat_1 = zeros(size(X_mat));
M_mat_2 = M_mat_1;
for idx = 1:freqpts
    M_mat_1(:,:,idx) = (Y_mat(:,:,idx)/X_mat(:,:,idx) + X_mat(:,:,idx)/Y_mat(:,:,idx))...
        /(X_mat(:,:,idx)/Y_mat(:,:,idx) - Y_mat(:,:,idx)/X_mat(:,:,idx));
    M_mat_2(:,:,idx) = (X_mat(:,:,idx)/Y_mat(:,:,idx) - Y_mat(:,:,idx)/X_mat(:,:,idx))\eye(2*numLines);
end

%% Calculate Zc,gamma
D_M1 = T_mat_1(numLines+1:end, numLines+1:end,:);
C_M2 = T_mat_2(numLines+1:end, 1:numLines,:);
gamma = acoshMat(D_M1)/LL; % Complex propagation constant((Nepers/m),(radians/m))
Zc = zeros(numLines,numLines,freqpts); % Characteristic line impedance(ohm)
for m = 1:freqpts
    sinh_gamma = funm(gamma(:,:,m)*LL,@sinh);
    if abs(sinh_gamma) < eps
        Zc(:,:,m) = eps;
    else
        Zc(:,:,m) = 0.5 * C_M2(:,:,m)\sinh_gamma;
    end
end
%% Calculate the RLGC matrices
rlgc_params = zGamma2rlgc(Zc,gamma,freq);

% One-line method
rlgc_params_oneLine = s2rlgc_mod2(SE8P_Data_1.S_Parameters,LL_1,freq,Z0);

%% Rebuilt S params
% Two-line method
s_params_rebuilt = snp2smp(sparameters(rlgc2s_mod(rlgc_params.R,rlgc_params.L,...
    rlgc_params.G,rlgc_params.C,LL_1,freq,Z0),freq,Z0),[1 2 3 4],Z0);

% One-line method
s_params_rebuilt_oneLine = snp2smp(sparameters(rlgc2s_mod(rlgc_params_oneLine.R,rlgc_params_oneLine.L,...
    rlgc_params_oneLine.G,rlgc_params_oneLine.C,LL_1,freq,Z0),freq,Z0),[1 2 3 4],Z0);

%% figure

figure('Name','s_params_rebuilt_oneLine');
rfplot(s_params_rebuilt_oneLine);
hold on
rfplot(s_params_1);
