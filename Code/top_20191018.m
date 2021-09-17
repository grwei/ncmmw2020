%% Basic information
% Given the 8-Port S parameters and port characteristic impedance of a MTL
% networks,
% extract the spice-Parameters:tabular RLGC model(W-element)
% Author Name: Wei Guorui
% Instructor Name: Xia Bin
% Created in: 2019-10-18 17:00
% Last modified: 2019-10-19 21:30
clc; clear; close all;
addpath(genpath('../code'))

%% Initialize
LL = 0.02; % the length of the transmission line in meters.
filename_4line = 'data/four_line/4line_4linetoRLCG_10201603.s8p';
SingleEnded8PortData = read(rfdata.data,filename_4line);
freq = SingleEnded8PortData.Freq;
numOfFreq = length(freq);
Z0 = SingleEnded8PortData.Z0;
ZL = SingleEnded8PortData.ZL;
ZS = SingleEnded8PortData.ZS;
SingleEnded8PortData.S_Parameters = snp2smp(SingleEnded8PortData.S_Parameters,...
    Z0,[4 8 2 6 1 5 3 7]);
numOfLines = size(SingleEnded8PortData.S_Parameters,1)/2;
S_params = sparameters(SingleEnded8PortData);


%% Compute the RLCG matrix
% 	           PORT NUMBERING (Modern style)
% 	    PORT#                               PORT#
% 	    Z1 = Z3 = 50.00                     Z2 = Z4 = 50.00
% 	           --------------------------
% 	          |                          |
% 	      1 --|                          |-- 2
% 	          |                          |
% 	     ...--|                          |-- ...
% 	          |                          |
% 	   2k-1 --|                          |-- 2k
% 	          |                          |
% 	   2N-1 --|                          |-- 2N
% 	          |                          |
% 	           --------------------------


% 1. Using Matlab RF Toolbox
[~, idx_nonpassive]= ispassive(S_params);
sparams_passive = makepassive(S_params);
% validate
[passivity, ~]= ispassive(sparams_passive);
err_S = max(abs((S_params.Parameters - sparams_passive.Parameters)./...
    (S_params.Parameters + eps)),[],[1 2 3]); % max passivity error
rlgc_params_matlab = s2rlgc_mod2(snp2smp(sparams_passive.Parameters,...
    Z0,[1 3 5 7 2 4 6 8]),LL,freq,Z0); % PORT NUMBERING (Classic style)

%% 2. Using New Method
% 2.1 Convert 2N-port S-parameters to (NC2) MMS Matrices
% Sxx_ij: mixed-mode S matrix of the ith and jth conductor
[smm_12.dd,smm_12.dc,smm_12.cd,smm_12.cc] = ...
    s2smm(snp2smp(S_params.Parameters,Z0,[1 2 3 4],Z0),1);
[smm_13.dd,smm_13.dc,smm_13.cd,smm_13.cc] = ...
    s2smm(snp2smp(S_params.Parameters,Z0,[1 2 5 6],Z0),1);
[smm_14.dd,smm_14.dc,smm_14.cd,smm_14.cc] = ...
    s2smm(snp2smp(S_params.Parameters,Z0,[1 2 7 8],Z0),1);
[smm_23.dd,smm_23.dc,smm_23.cd,smm_23.cc] = ...
    s2smm(snp2smp(S_params.Parameters,Z0,[3 4 5 6],Z0),1);
[smm_24.dd,smm_24.dc,smm_24.cd,smm_24.cc] = ...
    s2smm(snp2smp(S_params.Parameters,Z0,[3 4 7 8],Z0),1);
[smm_34.dd,smm_34.dc,smm_34.cd,smm_34.cc] = ...
    s2smm(snp2smp(S_params.Parameters,Z0,[5 6 7 8],Z0),1);

% 2.2 Calculate RLGC of each pairs
rlgc_12 = mms2rlgc(smm_12,LL,freq,Z0);
rlgc_13 = mms2rlgc(smm_13,LL,freq,Z0);
rlgc_14 = mms2rlgc(smm_14,LL,freq,Z0);
rlgc_23 = mms2rlgc(smm_23,LL,freq,Z0);
rlgc_24 = mms2rlgc(smm_24,LL,freq,Z0);
rlgc_34 = mms2rlgc(smm_34,LL,freq,Z0);

%% Verify each pairs
% Each pairs: Rebuilt S using RLGC extracted by RF Toolbox
s_params_rebuilt_matlab_12 = snp2smp(sparameters(rlgc2s_mod(...
    rlgc_params_matlab.R([1 2],[1 2],:),rlgc_params_matlab.L([1 2],[1 2],:),...
    rlgc_params_matlab.G([1 2],[1 2],:),rlgc_params_matlab.C([1 2],[1 2],:),...
    LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);
s_params_rebuilt_matlab_13 = snp2smp(sparameters(rlgc2s_mod(...
    rlgc_params_matlab.R([1 3],[1 3],:),rlgc_params_matlab.L([1 3],[1 3],:),...
    rlgc_params_matlab.G([1 3],[1 3],:),rlgc_params_matlab.C([1 3],[1 3],:),...
    LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);
s_params_rebuilt_matlab_14 = snp2smp(sparameters(rlgc2s_mod(...
    rlgc_params_matlab.R([1 4],[1 4],:),rlgc_params_matlab.L([1 4],[1 4],:),...
    rlgc_params_matlab.G([1 4],[1 4],:),rlgc_params_matlab.C([1 4],[1 4],:),...
    LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);
s_params_rebuilt_matlab_23 = snp2smp(sparameters(rlgc2s_mod(...
    rlgc_params_matlab.R([2 3],[2 3],:),rlgc_params_matlab.L([2 3],[2 3],:),...
    rlgc_params_matlab.G([2 3],[2 3],:),rlgc_params_matlab.C([2 3],[2 3],:),...
    LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);
s_params_rebuilt_matlab_24 = snp2smp(sparameters(rlgc2s_mod(...
    rlgc_params_matlab.R([2 4],[2 4],:),rlgc_params_matlab.L([2 4],[2 4],:),...
    rlgc_params_matlab.G([2 4],[2 4],:),rlgc_params_matlab.C([2 4],[2 4],:),...
    LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);
s_params_rebuilt_matlab_34 = snp2smp(sparameters(rlgc2s_mod(...
    rlgc_params_matlab.R([3 4],[3 4],:),rlgc_params_matlab.L([3 4],[3 4],:),...
    rlgc_params_matlab.G([3 4],[3 4],:),rlgc_params_matlab.C([3 4],[3 4],:),...
    LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);

% Each pairs: Rebuilt S using RLGC extracted by my method
s_params_rebuilt_12 = snp2smp(sparameters(rlgc2s_mod(rlgc_12.R,rlgc_12.L,...
    rlgc_12.G,rlgc_12.C,LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);
s_params_rebuilt_13 = snp2smp(sparameters(rlgc2s_mod(rlgc_13.R,rlgc_13.L,...
    rlgc_13.G,rlgc_13.C,LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);
s_params_rebuilt_14 = snp2smp(sparameters(rlgc2s_mod(rlgc_14.R,rlgc_14.L,...
    rlgc_14.G,rlgc_14.C,LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);
s_params_rebuilt_23 = snp2smp(sparameters(rlgc2s_mod(rlgc_23.R,rlgc_23.L,...
    rlgc_23.G,rlgc_23.C,LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);
s_params_rebuilt_24 = snp2smp(sparameters(rlgc2s_mod(rlgc_24.R,rlgc_24.L,...
    rlgc_24.G,rlgc_24.C,LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);
s_params_rebuilt_34 = snp2smp(sparameters(rlgc2s_mod(rlgc_34.R,rlgc_34.L,...
    rlgc_34.G,rlgc_34.C,LL,freq,Z0),freq,Z0),[1 3 2 4],Z0);

% 2.3 Combine RLGC of each pairs to general RLGC matrix
rlgc_params.R = zeros(numOfLines,numOfLines,numOfFreq);
rlgc_params.L = zeros(numOfLines,numOfLines,numOfFreq);
rlgc_params.G = zeros(numOfLines,numOfLines,numOfFreq);
rlgc_params.C = zeros(numOfLines,numOfLines,numOfFreq);
%% R matrix
% rlgc_params.R(1,1,:) =  rlgc_12.R(1,1,:) + rlgc_13.R(1,1,:) + rlgc_14.R(1,1,:);
% rlgc_params.R(2,2,:) =  rlgc_12.R(2,2,:) + rlgc_23.R(1,1,:) + rlgc_24.R(1,1,:);
% rlgc_params.R(3,3,:) =  rlgc_34.R(1,1,:) + rlgc_13.R(2,2,:) + rlgc_23.R(2,2,:);
% rlgc_params.R(4,4,:) =  rlgc_14.R(2,2,:) + rlgc_24.R(2,2,:) + rlgc_34.R(2,2,:);
% Debug
rlgc_params.R(1,1,:) =  rlgc_12.R(1,1,:);
rlgc_params.R(2,2,:) =  rlgc_12.R(2,2,:);
rlgc_params.R(3,3,:) =  rlgc_34.R(1,1,:);
rlgc_params.R(4,4,:) =  rlgc_34.R(2,2,:);

rlgc_params.R(1,2,:) =  rlgc_12.R(1,2,:);
rlgc_params.R(2,1,:) =  rlgc_12.R(2,1,:);
rlgc_params.R(3,4,:) =  rlgc_34.R(1,2,:);
rlgc_params.R(4,3,:) =  rlgc_34.R(2,1,:);

rlgc_params.R(1,3,:) =  rlgc_13.R(1,2,:);
rlgc_params.R(3,1,:) =  rlgc_13.R(2,1,:);
rlgc_params.R(2,4,:) =  rlgc_24.R(1,2,:);
rlgc_params.R(4,2,:) =  rlgc_24.R(2,1,:);
rlgc_params.R(1,4,:) =  rlgc_14.R(1,2,:);
rlgc_params.R(4,1,:) =  rlgc_14.R(2,1,:);
rlgc_params.R(2,3,:) =  rlgc_23.R(1,2,:);
rlgc_params.R(3,2,:) =  rlgc_23.R(2,1,:);
% Debug:
rlgc_params.R(1,3,:) =  eps;
rlgc_params.R(3,1,:) =  eps;
rlgc_params.R(2,4,:) =  eps;
rlgc_params.R(4,2,:) =  eps;
rlgc_params.R(1,4,:) =  eps;
rlgc_params.R(4,1,:) =  eps;
rlgc_params.R(2,3,:) =  eps;
rlgc_params.R(3,2,:) =  eps;
%% L matrix
% rlgc_params.L(1,1,:) =  rlgc_12.L(1,1,:) + rlgc_13.L(1,1,:) + rlgc_14.L(1,1,:);
% rlgc_params.L(2,2,:) =  rlgc_12.L(2,2,:) + rlgc_23.L(1,1,:) + rlgc_24.L(1,1,:);
% rlgc_params.L(3,3,:) =  rlgc_34.L(1,1,:) + rlgc_13.L(2,2,:) + rlgc_23.L(2,2,:);
% rlgc_params.L(4,4,:) =  rlgc_14.L(2,2,:) + rlgc_24.L(2,2,:) + rlgc_34.L(2,2,:);
% Debug:
rlgc_params.L(1,1,:) =  rlgc_12.L(1,1,:);
rlgc_params.L(2,2,:) =  rlgc_12.L(2,2,:);
rlgc_params.L(3,3,:) =  rlgc_34.L(1,1,:);
rlgc_params.L(4,4,:) =  rlgc_34.L(2,2,:);

rlgc_params.L(1,2,:) =  rlgc_12.L(1,2,:);
rlgc_params.L(2,1,:) =  rlgc_12.L(2,1,:);
rlgc_params.L(3,4,:) =  rlgc_34.L(1,2,:);
rlgc_params.L(4,3,:) =  rlgc_34.L(2,1,:);

rlgc_params.L(1,3,:) =  rlgc_13.L(1,2,:);
rlgc_params.L(3,1,:) =  rlgc_13.L(2,1,:);
rlgc_params.L(2,4,:) =  rlgc_24.L(1,2,:);
rlgc_params.L(4,2,:) =  rlgc_24.L(2,1,:);
rlgc_params.L(1,4,:) =  rlgc_14.L(1,2,:);
rlgc_params.L(4,1,:) =  rlgc_14.L(2,1,:);
rlgc_params.L(2,3,:) =  rlgc_23.L(1,2,:);
rlgc_params.L(3,2,:) =  rlgc_23.L(2,1,:);
% Debug:
rlgc_params.L(1,3,:) =  eps;
rlgc_params.L(3,1,:) =  eps;
rlgc_params.L(2,4,:) =  eps;
rlgc_params.L(4,2,:) =  eps;
rlgc_params.L(1,4,:) =  eps;
rlgc_params.L(4,1,:) =  eps;
rlgc_params.L(2,3,:) =  eps;
rlgc_params.L(3,2,:) =  eps;

%% C matrix
% rlgc_params.C(1,1,:) =  rlgc_12.C(1,1,:) + rlgc_13.C(1,1,:) + rlgc_14.C(1,1,:);
% rlgc_params.C(2,2,:) =  rlgc_12.C(2,2,:) + rlgc_23.C(1,1,:) + rlgc_24.C(1,1,:);
% rlgc_params.C(3,3,:) =  rlgc_34.C(1,1,:) + rlgc_13.C(2,2,:) + rlgc_23.C(2,2,:);
% rlgc_params.C(4,4,:) =  rlgc_14.C(2,2,:) + rlgc_24.C(2,2,:) + rlgc_34.C(2,2,:);
% Debug:
rlgc_params.C(1,1,:) =  rlgc_12.C(1,1,:);
rlgc_params.C(2,2,:) =  rlgc_12.C(2,2,:);
rlgc_params.C(3,3,:) =  rlgc_34.C(1,1,:);
rlgc_params.C(4,4,:) =  rlgc_14.C(2,2,:);

rlgc_params.C(1,2,:) =  rlgc_12.C(1,2,:);
rlgc_params.C(2,1,:) =  rlgc_12.C(2,1,:);
rlgc_params.C(3,4,:) =  rlgc_34.C(1,2,:);
rlgc_params.C(4,3,:) =  rlgc_34.C(2,1,:);

rlgc_params.C(1,3,:) =  rlgc_13.C(1,2,:);
rlgc_params.C(3,1,:) =  rlgc_13.C(2,1,:);
rlgc_params.C(2,4,:) =  rlgc_24.C(1,2,:);
rlgc_params.C(4,2,:) =  rlgc_24.C(2,1,:);
rlgc_params.C(1,4,:) =  rlgc_14.C(1,2,:);
rlgc_params.C(4,1,:) =  rlgc_14.C(2,1,:);
rlgc_params.C(2,3,:) =  rlgc_23.C(1,2,:);
rlgc_params.C(3,2,:) =  rlgc_23.C(2,1,:);
% Debug:
rlgc_params.C(1,3,:) =  eps;
rlgc_params.C(3,1,:) =  eps;
rlgc_params.C(2,4,:) =  eps;
rlgc_params.C(4,2,:) =  eps;
rlgc_params.C(1,4,:) =  eps;
rlgc_params.C(4,1,:) =  eps;
rlgc_params.C(2,3,:) =  eps;
rlgc_params.C(3,2,:) =  eps;
%% G matrix
% rlgc_params.G(1,1,:) =  rlgc_12.G(1,1,:) + rlgc_13.G(1,1,:) + rlgc_14.G(1,1,:);
% rlgc_params.G(2,2,:) =  rlgc_12.G(2,2,:) + rlgc_23.G(1,1,:) + rlgc_24.G(1,1,:);
% rlgc_params.G(3,3,:) =  rlgc_34.G(1,1,:) + rlgc_13.G(2,2,:) + rlgc_23.G(2,2,:);
% rlgc_params.G(4,4,:) =  rlgc_14.G(2,2,:) + rlgc_24.G(2,2,:) + rlgc_34.G(2,2,:);
% Debug:
rlgc_params.G(1,1,:) =  rlgc_12.G(1,1,:);
rlgc_params.G(2,2,:) =  rlgc_12.G(2,2,:);
rlgc_params.G(3,3,:) =  rlgc_34.G(1,1,:);
rlgc_params.G(4,4,:) =  rlgc_14.G(2,2,:);

rlgc_params.G(1,2,:) =  rlgc_12.G(1,2,:);
rlgc_params.G(2,1,:) =  rlgc_12.G(2,1,:);
rlgc_params.G(3,4,:) =  rlgc_34.G(1,2,:);
rlgc_params.G(4,3,:) =  rlgc_34.G(2,1,:);

rlgc_params.G(1,3,:) =  rlgc_13.G(1,2,:);
rlgc_params.G(3,1,:) =  rlgc_13.G(2,1,:);
rlgc_params.G(2,4,:) =  rlgc_24.G(1,2,:);
rlgc_params.G(4,2,:) =  rlgc_24.G(2,1,:);
rlgc_params.G(1,4,:) =  rlgc_14.G(1,2,:);
rlgc_params.G(4,1,:) =  rlgc_14.G(2,1,:);
rlgc_params.G(2,3,:) =  rlgc_23.G(1,2,:);
rlgc_params.G(3,2,:) =  rlgc_23.G(2,1,:);
% Debug:
rlgc_params.G(1,3,:) =  eps;
rlgc_params.G(3,1,:) =  eps;
rlgc_params.G(2,4,:) =  eps;
rlgc_params.G(4,2,:) =  eps;
rlgc_params.G(1,4,:) =  eps;
rlgc_params.G(4,1,:) =  eps;
rlgc_params.G(2,3,:) =  eps;
rlgc_params.G(3,2,:) =  eps;

%% Rebuilt S params
% RF Toolbox
s_params_rebuilt_matlab = snp2smp(sparameters(rlgc2s_mod(rlgc_params_matlab.R,rlgc_params_matlab.L,...
    rlgc_params_matlab.G,rlgc_params_matlab.C,LL,freq,Z0),freq,Z0),[1 5 2 6 3 7 4 8],Z0);
% My method
s_params_rebuilt = snp2smp(sparameters(rlgc2s_mod(rlgc_params.R,rlgc_params.L,...
    rlgc_params.G,rlgc_params.C,LL,freq,Z0),freq,Z0),[1 5 2 6 3 7 4 8],Z0);

%% figure

%% figure. Simulated Scaterring matrix
figure('Name','Simulated [S] Matrix')
subplot(2,4,1)
rfplot(S_params,1,1);
hold on
rfplot(S_params,2,2);
rfplot(S_params,7,7);
rfplot(S_params,8,8);
hold off
legend({'S11','S22','S77','S88'},'Location','best','NumColumns',2)
legend('boxoff')
title('Reflection: S11,S22,S77,S88')

subplot(2,4,2)
rfplot(S_params,3,3);
hold on
rfplot(S_params,4,4);
rfplot(S_params,5,5);
rfplot(S_params,6,6);
hold off
legend({'S33','S44','S55','S66'},'Location','best','NumColumns',2)
legend('boxoff')
title('Reflection: S33,S44,S55,S66')

subplot(2,4,3)
rfplot(S_params,1,2);
hold on
rfplot(S_params,2,1);
rfplot(S_params,7,8);
rfplot(S_params,8,7);
hold off
legend({'S12','S21','S78','S87'},'Location','best','NumColumns',2)
legend('boxoff')
title('Insertion: S12,S21,S78,S87')

subplot(2,4,4)
rfplot(S_params,3,4);
hold on
rfplot(S_params,4,3);
rfplot(S_params,5,6);
rfplot(S_params,6,5);
hold off
legend({'S34','S43','S56','S65'},'Location','best','NumColumns',2)
legend('boxoff')
title('Insertion: S34,S43,S56,S65')

subplot(2,4,5)
rfplot(S_params,1,4);
hold on
rfplot(S_params,4,1);
rfplot(S_params,2,3);
rfplot(S_params,3,2);
hold off
legend({'S14','S41','S23','S32'},'Location','best','NumColumns',2)
legend('boxoff')
title('Coupling: S14,S41,S23,S32')

subplot(2,4,6)
rfplot(S_params,5,8);
hold on
rfplot(S_params,8,5);
rfplot(S_params,7,6);
rfplot(S_params,6,7);
hold off
legend({'S58','S85','S76','S67'},'Location','best','NumColumns',2)
legend('boxoff')
title('Coupling: S58,S85,S76,S67')

subplot(2,4,7)
rfplot(S_params,1,3);
hold on
rfplot(S_params,3,1);
rfplot(S_params,2,4);
rfplot(S_params,4,2);
hold off
legend({'S13','S31','S24','S42'},'Location','best','NumColumns',2)
legend('boxoff')
title('Isolation: S13,S31,S24,S42')

subplot(2,4,8)
rfplot(S_params,5,7);
hold on
rfplot(S_params,7,5);
rfplot(S_params,6,8);
rfplot(S_params,8,6);
hold off
legend({'S57','S75','S68','S86'},'Location','best','NumColumns',2)
legend('boxoff')
title('Isolation: S57,S75,S68,S86')

% figure. S params assumes to zero
figure('Name','Simulated [S] Matrix: Assumes to zero')
subplot(232)
rfplot(S_params,1,6);
hold on
rfplot(S_params,6,1);
rfplot(S_params,5,2);
rfplot(S_params,2,5);
rfplot(S_params,3,8);
rfplot(S_params,8,3);
rfplot(S_params,7,4);
rfplot(S_params,4,7);
hold off
legend({'S16','S61','S52','S25','S38','S83','S74','S47'},'Location','best','NumColumns',3)
legend('boxoff')
title("Mid Neighbor: S16,S61,S52,S25," + newline + "S38,S83,S74,S47")

subplot(235)
rfplot(S_params,1,5);
hold on
rfplot(S_params,5,1);
rfplot(S_params,3,7);
rfplot(S_params,7,3);
rfplot(S_params,2,6);
rfplot(S_params,6,2);
rfplot(S_params,4,8);
rfplot(S_params,8,4);
hold off
legend({'S15','S51','S37','S73','S26','S62','S48','S84'},'Location','best','NumColumns',3)
legend('boxoff')
title("Mid Neighbor: S15,S51,S37,S73," + newline + "S26,S62,S48,S84")

subplot(233)
rfplot(S_params,1,8);
hold on
rfplot(S_params,8,1);
rfplot(S_params,2,7);
rfplot(S_params,7,2);
hold off
legend({'S18','S81','S27','S72'},'Location','best','NumColumns',2)
legend('boxoff')
title("Far Neighbor: S18,S81,S27,S72")

subplot(236)
rfplot(S_params,1,7);
hold on
rfplot(S_params,7,1);
rfplot(S_params,2,8);
rfplot(S_params,8,2);
hold off
legend({'S17','S71','S28','S82'},'Location','best','NumColumns',2)
legend('boxoff')
title("Far Neighbor: S17,S71,S28,S82")

subplot(231)
rfplot(S_params,3,6);
hold on
rfplot(S_params,6,3);
rfplot(S_params,4,5);
rfplot(S_params,5,4);
hold off
legend({'S36','S63','S45','S54'},'Location','best','NumColumns',2)
legend('boxoff')
title('Close Neighbor: S36,S63,S45,S54')

subplot(234)
rfplot(S_params,3,5);
hold on
rfplot(S_params,5,3);
rfplot(S_params,4,6);
rfplot(S_params,6,4);
hold off
legend({'S35','S53','S46','S64'},'Location','best','NumColumns',2)
legend('boxoff')
title('Close Neighbor: S35,S53,S46,S64')

%% figure. RLCG

%% matlab RF Toolbox

%% Validate the consistence of calculated RLGC of each pairs
figure('Name','RLGC: pairs Line 1j')
subplot(2,4,1)
plot(freq,squeeze(rlgc_12.R(1,1,:)))
hold on
plot(freq,squeeze(rlgc_13.R(1,1,:)))
plot(freq,squeeze(rlgc_14.R(1,1,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R11(Ohm/m)');
title('R11:Line 1j');
legend({'Line 12','Line 13','Line 14'},'Location','best')
legend('boxoff')

subplot(2,4,2)
plot(freq,squeeze(rlgc_12.L(1,1,:)))
hold on
plot(freq,squeeze(rlgc_13.L(1,1,:)))
plot(freq,squeeze(rlgc_14.L(1,1,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L11(H/m)');
title('L11:Line 1j');
legend({'Line 12','Line 13','Line 14'},'Location','best')
legend('boxoff')

subplot(2,4,3)
plot(freq,squeeze(rlgc_12.C(1,1,:)))
hold on
plot(freq,squeeze(rlgc_13.C(1,1,:)))
plot(freq,squeeze(rlgc_14.C(1,1,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C11(F/m)');
title('C11:Line 1j');
legend({'Line 12','Line 13','Line 14'},'Location','best')
legend('boxoff')

subplot(2,4,4)
plot(freq,squeeze(rlgc_12.G(1,1,:)))
hold on
plot(freq,squeeze(rlgc_13.G(1,1,:)))
plot(freq,squeeze(rlgc_14.G(1,1,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G22(S/m)');
title('G22:Line 1j');
legend({'Line 12','Line 13','Line 14'},'Location','best')
legend('boxoff')

subplot(2,4,5)
plot(freq,squeeze(rlgc_12.R(1,2,:)))
hold on
plot(freq,squeeze(rlgc_13.R(1,2,:)))
plot(freq,squeeze(rlgc_14.R(1,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R12(Ohm/m)');
title('R12:Line 1j');
legend({'Line 12','Line 13','Line 14'},'Location','best')
legend('boxoff')

subplot(2,4,6)
plot(freq,squeeze(rlgc_12.L(1,2,:)))
hold on
plot(freq,squeeze(rlgc_13.L(1,2,:)))
plot(freq,squeeze(rlgc_14.L(1,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L12(H/m)');
title('L12:Line 1j');
legend({'Line 12','Line 13','Line 14'},'Location','best')
legend('boxoff')

subplot(2,4,7)
plot(freq,squeeze(rlgc_12.C(1,2,:)))
hold on
plot(freq,squeeze(rlgc_13.C(1,2,:)))
plot(freq,squeeze(rlgc_14.C(1,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C12(F/m)');
title('C12:Line 1j');
legend({'Line 12','Line 13','Line 14'},'Location','best')
legend('boxoff')

subplot(2,4,8)
plot(freq,squeeze(rlgc_12.G(1,2,:)))
hold on
plot(freq,squeeze(rlgc_13.G(1,2,:)))
plot(freq,squeeze(rlgc_14.G(1,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G12(S/m)');
title('G12:Line 1j');
legend({'Line 12','Line 13','Line 14'},'Location','best')
legend('boxoff')

figure('Name','RLGC: pairs Line 2j')
subplot(2,4,1)
plot(freq,squeeze(rlgc_23.R(1,1,:)))
hold on
plot(freq,squeeze(rlgc_24.R(1,1,:)))
plot(freq,squeeze(rlgc_12.R(2,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R11(Ohm/m)');
title('R11:Line 2j');
legend({'Line 23','Line 24','Line 21'},'Location','best')
legend('boxoff')

subplot(2,4,2)
plot(freq,squeeze(rlgc_23.R(1,2,:)))
hold on
plot(freq,squeeze(rlgc_24.R(1,2,:)))
plot(freq,squeeze(rlgc_12.R(2,1,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R12(Ohm/m)');
title('R12:Line 2j');
legend({'Line 23','Line 24','Line 21'},'Location','best')
legend('boxoff')

subplot(2,4,3)
plot(freq,squeeze(rlgc_23.L(1,1,:)))
hold on
plot(freq,squeeze(rlgc_24.L(1,1,:)))
plot(freq,squeeze(rlgc_12.L(2,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L11(H/m)');
title('L11:Line 1j');
legend({'Line 23','Line 24','Line 21'},'Location','best')
legend('boxoff')

subplot(2,4,4)
plot(freq,squeeze(rlgc_23.L(1,2,:)))
hold on
plot(freq,squeeze(rlgc_24.L(1,2,:)))
plot(freq,squeeze(rlgc_12.L(2,1,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L12(Ohm/m)');
title('L12:Line 2j');
legend({'Line 23','Line 24','Line 21'},'Location','best')
legend('boxoff')

subplot(2,4,5)
plot(freq,squeeze(rlgc_23.C(1,1,:)))
hold on
plot(freq,squeeze(rlgc_24.C(1,1,:)))
plot(freq,squeeze(rlgc_12.C(2,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C11(F/m)');
title('C11:Line 2j');
legend({'Line 23','Line 24','Line 21'},'Location','best')
legend('boxoff')

subplot(2,4,6)
plot(freq,squeeze(rlgc_23.C(1,2,:)))
hold on
plot(freq,squeeze(rlgc_24.C(1,2,:)))
plot(freq,squeeze(rlgc_12.C(2,1,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C12(F/m)');
title('C12:Line 2j');
legend({'Line 23','Line 24','Line 21'},'Location','best')
legend('boxoff')

subplot(2,4,7)
plot(freq,squeeze(rlgc_23.G(1,1,:)))
hold on
plot(freq,squeeze(rlgc_24.G(1,1,:)))
plot(freq,squeeze(rlgc_12.G(2,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G11(S/m)');
title('G11:Line 2j');
legend({'Line 23','Line 24','Line 21'},'Location','best');
legend('boxoff')

subplot(2,4,8)
plot(freq,squeeze(rlgc_23.G(1,2,:)))
hold on
plot(freq,squeeze(rlgc_24.G(1,2,:)))
plot(freq,squeeze(rlgc_12.G(2,1,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G12(S/m)');
title('G12:Line 2j');
legend({'Line 23','Line 24','Line 21'},'Location','best');
legend('boxoff')

%% Matlab RF Toolbox "s2rlgc"
figure('Name','R Matrix: RF Toolbox')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params_matlab.R(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.R(2,2,:)))
plot(freq,squeeze(rlgc_params_matlab.R(3,3,:)))
plot(freq,squeeze(rlgc_params_matlab.R(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R11,R22,R33,R44');
legend({'R11','R22','R33','R44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params_matlab.R(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.R(2,1,:)))
plot(freq,squeeze(rlgc_params_matlab.R(3,4,:)))
plot(freq,squeeze(rlgc_params_matlab.R(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R12,R21,R34,R43');
legend({'R12','R21','R34','R43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params_matlab.R(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.R(3,1,:)))
plot(freq,squeeze(rlgc_params_matlab.R(2,4,:)))
plot(freq,squeeze(rlgc_params_matlab.R(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R13,R31,R24,R42');
legend({'R13','R31','R24','R42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params_matlab.R(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.R(4,1,:)))
plot(freq,squeeze(rlgc_params_matlab.R(2,3,:)))
plot(freq,squeeze(rlgc_params_matlab.R(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R14,R41,R23,R32');
legend({'R14','R41','R23','R32'},'Location','best','NumColumns',2)
legend('boxoff')

figure('Name','L Matrix: RF Toolbox')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params_matlab.L(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.L(2,2,:)))
plot(freq,squeeze(rlgc_params_matlab.L(3,3,:)))
plot(freq,squeeze(rlgc_params_matlab.L(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L11,L22,L33,L44');
legend({'L11','L22','L33','L44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params_matlab.L(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.L(2,1,:)))
plot(freq,squeeze(rlgc_params_matlab.L(3,4,:)))
plot(freq,squeeze(rlgc_params_matlab.L(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L12,L21,L34,L43');
legend({'L12','L21','L34','L43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params_matlab.L(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.L(3,1,:)))
plot(freq,squeeze(rlgc_params_matlab.L(2,4,:)))
plot(freq,squeeze(rlgc_params_matlab.L(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L13,L31,L24,L42');
legend({'L13','L31','L24','L42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params_matlab.L(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.L(4,1,:)))
plot(freq,squeeze(rlgc_params_matlab.L(2,3,:)))
plot(freq,squeeze(rlgc_params_matlab.L(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L14,L41,L23,L32');
legend({'L14','L41','L23','L32'},'Location','best','NumColumns',2)
legend('boxoff')

figure('Name','C Matrix: RF Toolbox')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params_matlab.C(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.C(2,2,:)))
plot(freq,squeeze(rlgc_params_matlab.C(3,3,:)))
plot(freq,squeeze(rlgc_params_matlab.C(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C11,C22,C33,C44');
legend({'C11','C22','C33','C44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params_matlab.C(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.C(2,1,:)))
plot(freq,squeeze(rlgc_params_matlab.C(3,4,:)))
plot(freq,squeeze(rlgc_params_matlab.C(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C12,C21,C34,C43');
legend({'C12','C21','C34','C43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params_matlab.C(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.C(3,1,:)))
plot(freq,squeeze(rlgc_params_matlab.C(2,4,:)))
plot(freq,squeeze(rlgc_params_matlab.C(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C13,C31,C24,C42');
legend({'C13','C31','C24','C42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params_matlab.C(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.C(4,1,:)))
plot(freq,squeeze(rlgc_params_matlab.C(2,3,:)))
plot(freq,squeeze(rlgc_params_matlab.C(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C14,C41,C23,C32');
legend({'C14','C41','C23','C32'},'Location','best','NumColumns',2)
legend('boxoff')

figure('Name','G Matrix: RF Toolbox')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params_matlab.G(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.G(2,2,:)))
plot(freq,squeeze(rlgc_params_matlab.G(3,3,:)))
plot(freq,squeeze(rlgc_params_matlab.G(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G11,G22,G33,G44');
legend({'G11','G22','G33','G44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params_matlab.G(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.G(2,1,:)))
plot(freq,squeeze(rlgc_params_matlab.G(3,4,:)))
plot(freq,squeeze(rlgc_params_matlab.G(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G12,G21,G34,G43');
legend({'G12','G21','G34','G43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params_matlab.G(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.G(3,1,:)))
plot(freq,squeeze(rlgc_params_matlab.G(2,4,:)))
plot(freq,squeeze(rlgc_params_matlab.G(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G13,G31,G24,G42');
legend({'G13','G31','G24','G42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params_matlab.G(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params_matlab.G(4,1,:)))
plot(freq,squeeze(rlgc_params_matlab.G(2,3,:)))
plot(freq,squeeze(rlgc_params_matlab.G(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G14,G41,G23,G32');
legend({'G14','G41','G23','G32'},'Location','best','NumColumns',2)
legend('boxoff')




%% My method
figure('Name','R Matrix: My method')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params.R(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params.R(2,2,:)))
plot(freq,squeeze(rlgc_params.R(3,3,:)))
plot(freq,squeeze(rlgc_params.R(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R11,R22,R33,R44');
legend({'R11','R22','R33','R44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params.R(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params.R(2,1,:)))
plot(freq,squeeze(rlgc_params.R(3,4,:)))
plot(freq,squeeze(rlgc_params.R(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R12,R21,R34,R43');
legend({'R12','R21','R34','R43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params.R(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params.R(3,1,:)))
plot(freq,squeeze(rlgc_params.R(2,4,:)))
plot(freq,squeeze(rlgc_params.R(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R13,R31,R24,R42');
legend({'R13','R31','R24','R42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params.R(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params.R(4,1,:)))
plot(freq,squeeze(rlgc_params.R(2,3,:)))
plot(freq,squeeze(rlgc_params.R(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R14,R41,R23,R32');
legend({'R14','R41','R23','R32'},'Location','best','NumColumns',2)
legend('boxoff')

figure('Name','L Matrix: My method')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params.L(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params.L(2,2,:)))
plot(freq,squeeze(rlgc_params.L(3,3,:)))
plot(freq,squeeze(rlgc_params.L(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L11,L22,L33,L44');
legend({'L11','L22','L33','L44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params.L(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params.L(2,1,:)))
plot(freq,squeeze(rlgc_params.L(3,4,:)))
plot(freq,squeeze(rlgc_params.L(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L12,L21,L34,L43');
legend({'L12','L21','L34','L43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params.L(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params.L(3,1,:)))
plot(freq,squeeze(rlgc_params.L(2,4,:)))
plot(freq,squeeze(rlgc_params.L(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L13,L31,L24,L42');
legend({'L13','L31','L24','L42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params.L(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params.L(4,1,:)))
plot(freq,squeeze(rlgc_params.L(2,3,:)))
plot(freq,squeeze(rlgc_params.L(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L14,L41,L23,L32');
legend({'L14','L41','L23','L32'},'Location','best','NumColumns',2)
legend('boxoff')

figure('Name','C Matrix: My method')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params.C(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params.C(2,2,:)))
plot(freq,squeeze(rlgc_params.C(3,3,:)))
plot(freq,squeeze(rlgc_params.C(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C11,C22,C33,C44');
legend({'C11','C22','C33','C44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params.C(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params.C(2,1,:)))
plot(freq,squeeze(rlgc_params.C(3,4,:)))
plot(freq,squeeze(rlgc_params.C(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C12,C21,C34,C43');
legend({'C12','C21','C34','C43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params.C(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params.C(3,1,:)))
plot(freq,squeeze(rlgc_params.C(2,4,:)))
plot(freq,squeeze(rlgc_params.C(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C13,C31,C24,C42');
legend({'C13','C31','C24','C42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params.C(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params.C(4,1,:)))
plot(freq,squeeze(rlgc_params.C(2,3,:)))
plot(freq,squeeze(rlgc_params.C(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C14,C41,C23,C32');
legend({'C14','C41','C23','C32'},'Location','best','NumColumns',2)
legend('boxoff')

figure('Name','G Matrix: My method')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params.G(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params.G(2,2,:)))
plot(freq,squeeze(rlgc_params.G(3,3,:)))
plot(freq,squeeze(rlgc_params.G(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G11,G22,G33,G44');
legend({'G11','G22','G33','G44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params.G(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params.G(2,1,:)))
plot(freq,squeeze(rlgc_params.G(3,4,:)))
plot(freq,squeeze(rlgc_params.G(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G12,G21,G34,G43');
legend({'G12','G21','G34','G43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params.G(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params.G(3,1,:)))
plot(freq,squeeze(rlgc_params.G(2,4,:)))
plot(freq,squeeze(rlgc_params.G(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G13,G31,G24,G42');
legend({'G13','G31','G24','G42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params.G(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params.G(4,1,:)))
plot(freq,squeeze(rlgc_params.G(2,3,:)))
plot(freq,squeeze(rlgc_params.G(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G14,G41,G23,G32');
legend({'G14','G41','G23','G32'},'Location','best','NumColumns',2)
legend('boxoff')

%% figure. Rebuilt Scaterring matrix

%% Using Matlab RF Toobox
figure('Name','Rebuilt [S] Matrix(Toolbox)')
subplot(2,4,1)
rfplot(s_params_rebuilt_matlab,1,1);
hold on
rfplot(s_params_rebuilt_matlab,2,2);
rfplot(s_params_rebuilt_matlab,7,7);
rfplot(s_params_rebuilt_matlab,8,8);
hold off
legend({'S11','S22','S77','S88'},'Location','best','NumColumns',2)
legend('boxoff')
title('Reflection: S11,S22,S77,S88')

subplot(2,4,2)
rfplot(s_params_rebuilt_matlab,3,3);
hold on
rfplot(s_params_rebuilt_matlab,4,4);
rfplot(s_params_rebuilt_matlab,5,5);
rfplot(s_params_rebuilt_matlab,6,6);
hold off
legend({'S33','S44','S55','S66'},'Location','best','NumColumns',2)
legend('boxoff')
title('Reflection: S33,S44,S55,S66')

subplot(2,4,3)
rfplot(s_params_rebuilt_matlab,1,2);
hold on
rfplot(s_params_rebuilt_matlab,2,1);
rfplot(s_params_rebuilt_matlab,7,8);
rfplot(s_params_rebuilt_matlab,8,7);
hold off
legend({'S12','S21','S78','S87'},'Location','best','NumColumns',2)
legend('boxoff')
title('Insertion: S12,S21,S78,S87')

subplot(2,4,4)
rfplot(s_params_rebuilt_matlab,3,4);
hold on
rfplot(s_params_rebuilt_matlab,4,3);
rfplot(s_params_rebuilt_matlab,5,6);
rfplot(s_params_rebuilt_matlab,6,5);
hold off
legend({'S34','S43','S56','S65'},'Location','best','NumColumns',2)
legend('boxoff')
title('Insertion: S34,S43,S56,S65')

subplot(2,4,5)
rfplot(s_params_rebuilt_matlab,1,4);
hold on
rfplot(s_params_rebuilt_matlab,4,1);
rfplot(s_params_rebuilt_matlab,2,3);
rfplot(s_params_rebuilt_matlab,3,2);
hold off
legend({'S14','S41','S23','S32'},'Location','best','NumColumns',2)
legend('boxoff')
title('Coupling: S14,S41,S23,S32')

subplot(2,4,6)
rfplot(s_params_rebuilt_matlab,5,8);
hold on
rfplot(s_params_rebuilt_matlab,8,5);
rfplot(s_params_rebuilt_matlab,7,6);
rfplot(s_params_rebuilt_matlab,6,7);
hold off
legend({'S58','S85','S76','S67'},'Location','best','NumColumns',2)
legend('boxoff')
title('Coupling: S58,S85,S76,S67')

subplot(2,4,7)
rfplot(s_params_rebuilt_matlab,1,3);
hold on
rfplot(s_params_rebuilt_matlab,3,1);
rfplot(s_params_rebuilt_matlab,2,4);
rfplot(s_params_rebuilt_matlab,4,2);
hold off
legend({'S13','S31','S24','S42'},'Location','best','NumColumns',2)
legend('boxoff')
title('Isolation: S13,S31,S24,S42')

subplot(2,4,8)
rfplot(s_params_rebuilt_matlab,5,7);
hold on
rfplot(s_params_rebuilt_matlab,7,5);
rfplot(s_params_rebuilt_matlab,6,8);
rfplot(s_params_rebuilt_matlab,8,6);
hold off
legend({'S57','S75','S68','S86'},'Location','best','NumColumns',2)
legend('boxoff')
title('Isolation: S57,S75,S68,S86')

% figure. S params assumes to zero
figure('Name','Rebuilt [S] Matrix(Toolbox): Assumes to zero')
subplot(232)
rfplot(s_params_rebuilt_matlab,1,6);
hold on
rfplot(s_params_rebuilt_matlab,6,1);
rfplot(s_params_rebuilt_matlab,5,2);
rfplot(s_params_rebuilt_matlab,2,5);
rfplot(s_params_rebuilt_matlab,3,8);
rfplot(s_params_rebuilt_matlab,8,3);
rfplot(s_params_rebuilt_matlab,7,4);
rfplot(s_params_rebuilt_matlab,4,7);
hold off
legend({'S16','S61','S52','S25','S38','S83','S74','S47'},'Location','best','NumColumns',2)
legend('boxoff')
title("Mid Neighbor: S16,S61,S52,S25," + newline + "S38,S83,S74,S47")

subplot(235)
rfplot(s_params_rebuilt_matlab,1,5);
hold on
rfplot(s_params_rebuilt_matlab,5,1);
rfplot(s_params_rebuilt_matlab,3,7);
rfplot(s_params_rebuilt_matlab,7,3);
rfplot(s_params_rebuilt_matlab,2,6);
rfplot(s_params_rebuilt_matlab,6,2);
rfplot(s_params_rebuilt_matlab,4,8);
rfplot(s_params_rebuilt_matlab,8,4);
hold off
legend({'S15','S51','S37','S73','S26','S62','S48','S84'},'Location','best','NumColumns',2)
legend('boxoff')
title("Mid Neighbor: S15,S51,S37,S73," + newline + "S26,S62,S48,S84")

subplot(233)
rfplot(s_params_rebuilt_matlab,1,8);
hold on
rfplot(s_params_rebuilt_matlab,8,1);
rfplot(s_params_rebuilt_matlab,2,7);
rfplot(s_params_rebuilt_matlab,7,2);
hold off
legend({'S18','S81','S27','S72'},'Location','best','NumColumns',2)
legend('boxoff')
title("Far Neighbor: S18,S81,S27,S72")

subplot(236)
rfplot(s_params_rebuilt_matlab,1,7);
hold on
rfplot(s_params_rebuilt_matlab,7,1);
rfplot(s_params_rebuilt_matlab,2,8);
rfplot(s_params_rebuilt_matlab,8,2);
hold off
legend({'S17','S71','S28','S82'},'Location','best','NumColumns',2)
legend('boxoff')
title("Far Neighbor: S17,S71,S28,S82")

subplot(231)
rfplot(s_params_rebuilt_matlab,3,6);
hold on
rfplot(s_params_rebuilt_matlab,6,3);
rfplot(s_params_rebuilt_matlab,4,5);
rfplot(s_params_rebuilt_matlab,5,4);
hold off
legend({'S36','S63','S45','S54'},'Location','best','NumColumns',2)
legend('boxoff')
title('Close Neighbor: S36,S63,S45,S54')

subplot(234)
rfplot(s_params_rebuilt_matlab,3,5);
hold on
rfplot(s_params_rebuilt_matlab,5,3);
rfplot(s_params_rebuilt_matlab,4,6);
rfplot(s_params_rebuilt_matlab,6,4);
hold off
legend({'S35','S53','S46','S64'},'Location','best','NumColumns',2)
legend('boxoff')
title('Close Neighbor: S35,S53,S46,S64')

%% Using My Method
figure('Name','Rebuilt [S] Matrix(My method)')
subplot(2,4,1)
rfplot(s_params_rebuilt,1,1);
hold on
rfplot(s_params_rebuilt,2,2);
rfplot(s_params_rebuilt,7,7);
rfplot(s_params_rebuilt,8,8);
hold off
legend({'S11','S22','S77','S88'},'Location','best','NumColumns',2)
legend('boxoff')
title('Reflection: S11,S22,S77,S88')

subplot(2,4,2)
rfplot(s_params_rebuilt,3,3);
hold on
rfplot(s_params_rebuilt,4,4);
rfplot(s_params_rebuilt,5,5);
rfplot(s_params_rebuilt,6,6);
hold off
legend({'S33','S44','S55','S66'},'Location','best','NumColumns',2)
legend('boxoff')
title('Reflection: S33,S44,S55,S66')

subplot(2,4,3)
rfplot(s_params_rebuilt,1,2);
hold on
rfplot(s_params_rebuilt,2,1);
rfplot(s_params_rebuilt,7,8);
rfplot(s_params_rebuilt,8,7);
hold off
legend({'S12','S21','S78','S87'},'Location','best','NumColumns',2)
legend('boxoff')
title('Insertion: S12,S21,S78,S87')

subplot(2,4,4)
rfplot(s_params_rebuilt,3,4);
hold on
rfplot(s_params_rebuilt,4,3);
rfplot(s_params_rebuilt,5,6);
rfplot(s_params_rebuilt,6,5);
hold off
legend({'S34','S43','S56','S65'},'Location','best','NumColumns',2)
legend('boxoff')
title('Insertion: S34,S43,S56,S65')

subplot(2,4,5)
rfplot(s_params_rebuilt,1,4);
hold on
rfplot(s_params_rebuilt,4,1);
rfplot(s_params_rebuilt,2,3);
rfplot(s_params_rebuilt,3,2);
hold off
legend({'S14','S41','S23','S32'},'Location','best','NumColumns',2)
legend('boxoff')
title('Coupling: S14,S41,S23,S32')

subplot(2,4,6)
rfplot(s_params_rebuilt,5,8);
hold on
rfplot(s_params_rebuilt,8,5);
rfplot(s_params_rebuilt,7,6);
rfplot(s_params_rebuilt,6,7);
hold off
legend({'S58','S85','S76','S67'},'Location','best','NumColumns',2)
legend('boxoff')
title('Coupling: S58,S85,S76,S67')

subplot(2,4,7)
rfplot(s_params_rebuilt,1,3);
hold on
rfplot(s_params_rebuilt,3,1);
rfplot(s_params_rebuilt,2,4);
rfplot(s_params_rebuilt,4,2);
hold off
legend({'S13','S31','S24','S42'},'Location','best','NumColumns',2)
legend('boxoff')
title('Isolation: S13,S31,S24,S42')

subplot(2,4,8)
rfplot(s_params_rebuilt,5,7);
hold on
rfplot(s_params_rebuilt,7,5);
rfplot(s_params_rebuilt,6,8);
rfplot(s_params_rebuilt,8,6);
hold off
legend({'S57','S75','S68','S86'},'Location','best','NumColumns',2)
legend('boxoff')
title('Isolation: S57,S75,S68,S86')

% figure. S params assumes to zero
figure('Name','Rebuilt [S] Matrix(Toolbox): Assumes to zero')
subplot(232)
rfplot(s_params_rebuilt,1,6);
hold on
rfplot(s_params_rebuilt,6,1);
rfplot(s_params_rebuilt,5,2);
rfplot(s_params_rebuilt,2,5);
rfplot(s_params_rebuilt,3,8);
rfplot(s_params_rebuilt,8,3);
rfplot(s_params_rebuilt,7,4);
rfplot(s_params_rebuilt,4,7);
hold off
legend({'S16','S61','S52','S25','S38','S83','S74','S47'},'Location','best','NumColumns',2)
legend('boxoff')
title("Mid Neighbor: S16,S61,S52,S25," + newline + "S38,S83,S74,S47")

subplot(235)
rfplot(s_params_rebuilt,1,5);
hold on
rfplot(s_params_rebuilt,5,1);
rfplot(s_params_rebuilt,3,7);
rfplot(s_params_rebuilt,7,3);
rfplot(s_params_rebuilt,2,6);
rfplot(s_params_rebuilt,6,2);
rfplot(s_params_rebuilt,4,8);
rfplot(s_params_rebuilt,8,4);
hold off
legend({'S15','S51','S37','S73','S26','S62','S48','S84'},'Location','best','NumColumns',2)
legend('boxoff')
title("Mid Neighbor: S15,S51,S37,S73," + newline + "S26,S62,S48,S84")

subplot(233)
rfplot(s_params_rebuilt,1,8);
hold on
rfplot(s_params_rebuilt,8,1);
rfplot(s_params_rebuilt,2,7);
rfplot(s_params_rebuilt,7,2);
hold off
legend({'S18','S81','S27','S72'},'Location','best','NumColumns',2)
legend('boxoff')
title("Far Neighbor: S18,S81,S27,S72")

subplot(236)
rfplot(s_params_rebuilt,1,7);
hold on
rfplot(s_params_rebuilt,7,1);
rfplot(s_params_rebuilt,2,8);
rfplot(s_params_rebuilt,8,2);
hold off
legend({'S17','S71','S28','S82'},'Location','best','NumColumns',2)
legend('boxoff')
title("Far Neighbor: S17,S71,S28,S82")

subplot(231)
rfplot(s_params_rebuilt,3,6);
hold on
rfplot(s_params_rebuilt,6,3);
rfplot(s_params_rebuilt,4,5);
rfplot(s_params_rebuilt,5,4);
hold off
legend({'S36','S63','S45','S54'},'Location','best','NumColumns',2)
legend('boxoff')
title('Close Neighbor: S36,S63,S45,S54')

subplot(234)
rfplot(s_params_rebuilt,3,5);
hold on
rfplot(s_params_rebuilt,5,3);
rfplot(s_params_rebuilt,4,6);
rfplot(s_params_rebuilt,6,4);
hold off
legend({'S35','S53','S46','S64'},'Location','best','NumColumns',2)
legend('boxoff')
title('Close Neighbor: S35,S53,S46,S64')

%% Verify each pairs RLCG

%% My method
figure('Name','Verify each pairs: s_params_rebuilt(My method)')
subplot(2,3,1)
rfplot(s_params_rebuilt_12,1,1); % S11 = S22 = S33 = S44
hold on
rfplot(s_params_rebuilt_12,1,2); % S12 = S21 = S34 = S43
rfplot(s_params_rebuilt_12,1,3); % S13 = S31 = S24 = S42
rfplot(s_params_rebuilt_12,1,4); % S14 = S41 = S23 = S32
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 12')

subplot(2,3,2)
rfplot(s_params_rebuilt_13,1,1);
hold on
rfplot(s_params_rebuilt_13,1,2);
rfplot(s_params_rebuilt_13,1,3);
rfplot(s_params_rebuilt_13,1,4);
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 13')

subplot(2,3,3)
rfplot(s_params_rebuilt_14,1,1);
hold on
rfplot(s_params_rebuilt_14,1,2);
rfplot(s_params_rebuilt_14,1,3);
rfplot(s_params_rebuilt_14,1,4);
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 14')

subplot(2,3,4)
rfplot(s_params_rebuilt_23,1,1);
hold on
rfplot(s_params_rebuilt_23,1,2);
rfplot(s_params_rebuilt_23,1,3);
rfplot(s_params_rebuilt_23,1,4);
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 23')

subplot(2,3,5)
rfplot(s_params_rebuilt_24,1,1);
hold on
rfplot(s_params_rebuilt_24,1,2);
rfplot(s_params_rebuilt_24,1,3);
rfplot(s_params_rebuilt_24,1,4);
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 24')

subplot(2,3,6)
rfplot(s_params_rebuilt_34,1,1);
hold on
rfplot(s_params_rebuilt_34,1,2);
rfplot(s_params_rebuilt_34,1,3);
rfplot(s_params_rebuilt_34,1,4);
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 34')

%% RF Toolbox
%% My method
figure('Name','Verify each pairs: s_params_rebuilt(Toolbox)')
subplot(2,3,1)
rfplot(s_params_rebuilt_matlab_12,1,1); % S11 = S22 = S33 = S44
hold on
rfplot(s_params_rebuilt_matlab_12,1,2); % S12 = S21 = S34 = S43
rfplot(s_params_rebuilt_matlab_12,1,3); % S13 = S31 = S24 = S42
rfplot(s_params_rebuilt_matlab_12,1,4); % S14 = S41 = S23 = S32
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 12')

subplot(2,3,2)
rfplot(s_params_rebuilt_matlab_13,1,1);
hold on
rfplot(s_params_rebuilt_matlab_13,1,2);
rfplot(s_params_rebuilt_matlab_13,1,3);
rfplot(s_params_rebuilt_matlab_13,1,4);
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 13')

subplot(2,3,3)
rfplot(s_params_rebuilt_matlab_14,1,1);
hold on
rfplot(s_params_rebuilt_matlab_14,1,2);
rfplot(s_params_rebuilt_matlab_14,1,3);
rfplot(s_params_rebuilt_matlab_14,1,4);
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 14')

subplot(2,3,4)
rfplot(s_params_rebuilt_matlab_23,1,1);
hold on
rfplot(s_params_rebuilt_matlab_23,1,2);
rfplot(s_params_rebuilt_matlab_23,1,3);
rfplot(s_params_rebuilt_matlab_23,1,4);
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 23')

subplot(2,3,5)
rfplot(s_params_rebuilt_matlab_24,1,1);
hold on
rfplot(s_params_rebuilt_matlab_24,1,2);
rfplot(s_params_rebuilt_matlab_24,1,3);
rfplot(s_params_rebuilt_matlab_24,1,4);
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 24')

subplot(2,3,6)
rfplot(s_params_rebuilt_matlab_34,1,1);
hold on
rfplot(s_params_rebuilt_matlab_34,1,2);
rfplot(s_params_rebuilt_matlab_34,1,3);
rfplot(s_params_rebuilt_matlab_34,1,4);
hold off
% legend('off')
legend({'S11','S12','S13','S14'},'Location','best','NumColumns',2)
legend('boxoff')
title('s params rebuilt 34')

%% Verify: is submatrix: RLGC -> S?

%% Line 12
% Toolbox
figure('Name','is submatrix? Line 12(Toolbox)')
subplot(2,2,1)
rfplot(s_params_rebuilt_matlab_12,1,1); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt_matlab,1,1); % submat of rebuilt S(from total rlgc)
rfplot(S_params,1,1) % submat of simulated S
hold off
% legend('off')
legend({'S11-rlgcSub','S11-rlgcTot','S11-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S11(Toolbox)')

subplot(2,2,2)
rfplot(s_params_rebuilt_matlab_12,1,2); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt_matlab,1,2); % submat of rebuilt S(from total rlgc)
rfplot(S_params,1,2) % submat of simulated S
hold off
% legend('off')
legend({'S12-rlgcSub','S12-rlgcTot','S12-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S12(Toolbox)')

subplot(2,2,3)
rfplot(s_params_rebuilt_matlab_12,1,3); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt_matlab,1,3); % submat of rebuilt S(from total rlgc)
rfplot(S_params,1,3) % submat of simulated S
hold off
% legend('off')
legend({'S13-rlgcSub','S13-rlgcTot','S13-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S13(Toolbox)')

subplot(2,2,4)
rfplot(s_params_rebuilt_matlab_12,1,4); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt_matlab,1,4); % submat of rebuilt S(from total rlgc)
rfplot(S_params,1,4) % submat of simulated S
hold off
% legend('off')
legend({'S14-rlgcSub','S14-rlgcTot','S14-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S14(Toolbox)')

% My method
figure('Name','is submatrix? Line 12(My method)')
subplot(2,2,1)
rfplot(s_params_rebuilt_12,1,1); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt,1,1); % submat of rebuilt S(from total rlgc)
rfplot(S_params,1,1) % submat of simulated S
hold off
% legend('off')
legend({'S11-rlgcSub','S11-rlgcTot','S11-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S11(My method)')

subplot(2,2,2)
rfplot(s_params_rebuilt_12,1,2); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt,1,2); % submat of rebuilt S(from total rlgc)
rfplot(S_params,1,2) % submat of simulated S
hold off
% legend('off')
legend({'S12-rlgcSub','S12-rlgcTot','S12-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S12(My method)')

subplot(2,2,3)
rfplot(s_params_rebuilt_12,1,3); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt,1,3); % submat of rebuilt S(from total rlgc)
rfplot(S_params,1,3) % submat of simulated S
hold off
% legend('off')
legend({'S13-rlgcSub','S13-rlgcTot','S13-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S13(My method)')

subplot(2,2,4)
rfplot(s_params_rebuilt_12,1,4); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt,1,4); % submat of rebuilt S(from total rlgc)
rfplot(S_params,1,4) % submat of simulated S
hold off
% legend('off')
legend({'S14-rlgcSub','S14-rlgcTot','S14-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S14(My method)')

%% Line 34
% Toolbox
figure('Name','is submatrix? Line 34(Toolbox)')
subplot(2,2,1)
rfplot(s_params_rebuilt_matlab_34,1,1); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt_matlab,5,5); % submat of rebuilt S(from total rlgc)
rfplot(S_params,5,5) % submat of simulated S
hold off
% legend('off')
legend({'S55-rlgcSub','S55-rlgcTot','S55-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S55(Toolbox)')

subplot(2,2,2)
rfplot(s_params_rebuilt_matlab_34,1,2); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt_matlab,5,6); % submat of rebuilt S(from total rlgc)
rfplot(S_params,5,6) % submat of simulated S
hold off
% legend('off')
legend({'S56-rlgcSub','S56-rlgcTot','S56-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S56(Toolbox)')

subplot(2,2,3)
rfplot(s_params_rebuilt_matlab_34,1,3); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt_matlab,5,7); % submat of rebuilt S(from total rlgc)
rfplot(S_params,5,7) % submat of simulated S
hold off
% legend('off')
legend({'S57-rlgcSub','S57-rlgcTot','S57-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S57(Toolbox)')

subplot(2,2,4)
rfplot(s_params_rebuilt_matlab_34,1,4); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt_matlab,5,8); % submat of rebuilt S(from total rlgc)
rfplot(S_params,5,8) % submat of simulated S
hold off
% legend('off')
legend({'S58-rlgcSub','S58-rlgcTot','S58-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S58(Toolbox)')

% My method
figure('Name','is submatrix? Line 34(My method)')
subplot(2,2,1)
rfplot(s_params_rebuilt_34,1,1); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt,5,5); % submat of rebuilt S(from total rlgc)
rfplot(S_params,5,5) % submat of simulated S
hold off
% legend('off')
legend({'S55-rlgcSub','S55-rlgcTot','S55-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S55(My method)')

subplot(2,2,2)
rfplot(s_params_rebuilt_34,1,2); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt,5,6); % submat of rebuilt S(from total rlgc)
rfplot(S_params,5,6) % submat of simulated S
hold off
% legend('off')
legend({'S56-rlgcSub','S56-rlgcTot','S56-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S56(My method)')

subplot(2,2,3)
rfplot(s_params_rebuilt_34,1,3); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt,5,7); % submat of rebuilt S(from total rlgc)
rfplot(S_params,5,7) % submat of simulated S
hold off
% legend('off')
legend({'S57-rlgcSub','S57-rlgcTot','S57-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S57(My method)')

subplot(2,2,4)
rfplot(s_params_rebuilt_34,1,4); % submat of RLGC -> S
hold on
rfplot(s_params_rebuilt,5,8); % submat of rebuilt S(from total rlgc)
rfplot(S_params,5,8) % submat of simulated S
hold off
% legend('off')
legend({'S58-rlgcSub','S58-rlgcTot','S58-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('is submatrix? S58(My method)')
