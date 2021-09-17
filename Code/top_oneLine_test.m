%% Basic information
% The verification of the following paper(Section II):
% M. K. Sampath, "On addressing the practical issues in the extraction of RLGC parameters for lossy multiconductor transmission lines using S-parameter models," 2008 IEEE-EPEP Electrical Performance of Electronic Packaging, San Jose, CA, 2008, pp. 259-262.
% Author Name: Wei Guorui
% Instructor Name: Xia Bin
% Created at: 2019-10-25 22:40
% Last modified: 2019-10-25 
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
linelength = 0.0224; % line length(meters).
filename_2line = 'data/SI/d4line_22_4_mm_201910242117.s4p';
SE4P_Data = read(rfdata.data,filename_2line); % SingleEnded8PortData
freq = SE4P_Data.Freq;
freqpts = length(freq);
Z0 = SE4P_Data.Z0;
ZL = SE4P_Data.ZL;
ZS = SE4P_Data.ZS;
SE4P_Data.S_Parameters = snp2smp(SE4P_Data.S_Parameters,...
    Z0,[1 2 3 4]);
numLines = size(SE4P_Data.S_Parameters,1)/2;
s_params = sparameters(SE4P_Data);
S11_simu = squeeze(SE4P_Data.S_Parameters(1,1,:));
%% Convert 2N¡Á2N S matrix to 2N¡Á2N T Matrix(ABCD) and Z matrix
T_mat = s2t_mod(SE4P_Data.S_Parameters, Z0); % notice the difference in "s2rlgc": S -> T.
Td = T_mat(numLines+1:end, numLines+1:end,:);
% Td(:,:,idx) = ((I-s11(:,:,idx))*(I+s22(:,:,idx))+s12(:,:,idx)*s21(:,:,idx))/(2*s21(:,:,idx));
Z_params = s2z(SE4P_Data.S_Parameters,Z0);

%% Calculate Zc,gamma
gamma = acoshMat(Td)/linelength; % Complex propagation constant((Nepers/m),(radians/m))
Zc = zeros(numLines,numLines,freqpts); % Characteristic line impedance(ohm)
for m = 1:freqpts
    sinh_gamma = funm(gamma(:,:,m)*linelength,@sinh);
    if abs(sinh_gamma) < eps
        Zc(:,:,m) = eps;
    else
        Zc(:,:,m) = Z_params(numLines+1:end,1:numLines,m)*sinh_gamma;
    end
end

%% Calculate the RLGC matrices
rlgc_params = zGamma2rlgc(Zc,gamma,freq);

% Compare with RF Toolbox("s2rlgc")
rlgc_params_matlab= s2rlgc_mod2(SE4P_Data.S_Parameters,linelength,freq,Z0);

%% Rebuilt S params
% One-line method
s_params_rebuilt = snp2smp(sparameters(rlgc2s_mod(rlgc_params.R,rlgc_params.L,...
    rlgc_params.G,rlgc_params.C,linelength,freq,Z0),freq,Z0),[1 2 3 4],Z0);

% Compare with RF Toolbox("s2rlgc")
s_params_rebuilt_matlab = snp2smp(sparameters(rlgc2s_mod(rlgc_params_matlab.R,rlgc_params_matlab.L,...
    rlgc_params_matlab.G,rlgc_params_matlab.C,linelength,freq,Z0),freq,Z0),[1 2 3 4],Z0);

%% Fitting the RLCG parameters obtained by the closed-form equations with the model
% outliers = find((abs(imag(Z_c)) > tol_Zc_imag) | (freq' > freq_inlier_high));
% [R_fitresult, R_gof] = R_fit(freq', R, outliers);

%% figure.

%% S-params simu