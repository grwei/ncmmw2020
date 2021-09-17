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
LL_1 = 0.020; % line length(meters).
LL_2 = 0.015; 
LL = LL_1 - LL_2;
filename_TL1 = 'data/four_line/d4Line_20mm_201910232329.s8p';
filename_TL2 = 'data/four_line/d4Line_15mm_201910232329.s8p';
SE8P_Data_1 = read(rfdata.data,filename_TL1); % SingleEnded8PortData
SE8P_Data_2 = read(rfdata.data,filename_TL2); % SingleEnded8PortData
Z0 = SE8P_Data_1.Z0;
ZL = SE8P_Data_1.ZL;
ZS = SE8P_Data_1.ZS;
SE8P_Data_1.S_Parameters = snp2smp(SE8P_Data_1.S_Parameters,...
    Z0,[4 2 1 3 8 6 5 7]);
SE8P_Data_2.S_Parameters = snp2smp(SE8P_Data_2.S_Parameters,...
    Z0,[4 2 1 3 8 6 5 7]);
freq = SE8P_Data_1.Freq;
freqpts = length(freq);
numLines = size(SE8P_Data_1.S_Parameters,1)/2;
s_params_1 = sparameters(SE8P_Data_1);
s_params_2 = sparameters(SE8P_Data_2);

%% Convert 2N¡Á2N S matrix to 2N¡Á2N T Matrix(ABCD) and Z matrix
T_mat_1 = s2abcd(SE8P_Data_1.S_Parameters, Z0); % notice the difference in "s2rlgc": S -> T.
T_mat_2 = s2abcd(SE8P_Data_2.S_Parameters, Z0); % notice the difference in "s2rlgc": S -> T.

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

%% Rebuilt S params
% Two-line method
s_params_rebuilt = snp2smp(sparameters(rlgc2s_mod(rlgc_params.R,rlgc_params.L,...
    rlgc_params.G,rlgc_params.C,LL_1,freq,Z0),freq,Z0),[1 2 3 4 5 6 7 8],Z0);

%% Figure.

%% Extracted RLGC params: Two-line method
figure('Name','R Matrix: Two-line method')
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

%% Verify Rebuilt S params(Two line method)

%% Line 12
figure('Name','Verify Rebuilt S params: Line 12')
subplot(2,2,1)
rfplot(s_params_rebuilt,1,1); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,1) % submat of simulated S
hold off
% legend('off')
legend({'S11-rebuilt','S11-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S11')

subplot(2,2,2)
rfplot(s_params_rebuilt,1,2); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,2) % submat of simulated S
hold off
% legend('off')
legend({'S12-rebuilt','S12-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S12')

subplot(2,2,3)
rfplot(s_params_rebuilt,1,5); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,5) % submat of simulated S
hold off
% legend('off')
legend({'S15-rebuilt','S15-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S15')

subplot(2,2,4)
rfplot(s_params_rebuilt,1,6); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,6) % submat of simulated S
hold off
% legend('off')
legend({'S16-rebuilt','S16-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S16')

%% Line 13
figure('Name','Verify Rebuilt S params: Line 13')
subplot(2,2,1)
rfplot(s_params_rebuilt,1,1); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,1) % submat of simulated S
hold off
% legend('off')
legend({'S11-rebuilt','S11-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S11')

subplot(2,2,2)
rfplot(s_params_rebuilt,1,3); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,3) % submat of simulated S
hold off
% legend('off')
legend({'S13-rebuilt','S13-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S13')

subplot(2,2,3)
rfplot(s_params_rebuilt,1,5); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,5) % submat of simulated S
hold off
% legend('off')
legend({'S15-rebuilt','S15-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S15')

subplot(2,2,4)
rfplot(s_params_rebuilt,1,7); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,7) % submat of simulated S
hold off
% legend('off')
legend({'S17-rebuilt','S17-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S17')

%% Line 14
figure('Name','Verify Rebuilt S params: Line 14')
subplot(2,2,1)
rfplot(s_params_rebuilt,1,1); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,1) % submat of simulated S
hold off
% legend('off')
legend({'S11-rebuilt','S11-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S11')

subplot(2,2,2)
rfplot(s_params_rebuilt,1,4); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,4) % submat of simulated S
hold off
% legend('off')
legend({'S14-rebuilt','S14-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S14')

subplot(2,2,3)
rfplot(s_params_rebuilt,1,5); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,5) % submat of simulated S
hold off
% legend('off')
legend({'S15-rebuilt','S15-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S15')

subplot(2,2,4)
rfplot(s_params_rebuilt,1,8); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,1,8) % submat of simulated S
hold off
% legend('off')
legend({'S18-rebuilt','S18-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S18')

%% Line 23
figure('Name','Verify Rebuilt S params: Line 23')
subplot(2,2,1)
rfplot(s_params_rebuilt,2,2); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,2,2) % submat of simulated S
hold off
% legend('off')
legend({'S22-rebuilt','S22-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S22')

subplot(2,2,2)
rfplot(s_params_rebuilt,2,6); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,2,6) % submat of simulated S
hold off
% legend('off')
legend({'S26-rebuilt','S26-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S26')

subplot(2,2,3)
rfplot(s_params_rebuilt,2,3); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,2,3) % submat of simulated S
hold off
% legend('off')
legend({'S23-rebuilt','S23-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S23')

subplot(2,2,4)
rfplot(s_params_rebuilt,2,7); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,2,7) % submat of simulated S
hold off
% legend('off')
legend({'S27-rebuilt','S27-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S27')

%% Line 24
figure('Name','Verify Rebuilt S params: Line 24')
subplot(2,2,1)
rfplot(s_params_rebuilt,2,2); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,2,2) % submat of simulated S
hold off
% legend('off')
legend({'S11-rebuilt','S11-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S22')

subplot(2,2,2)
rfplot(s_params_rebuilt,2,6); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,2,6) % submat of simulated S
hold off
% legend('off')
legend({'S26-rebuilt','S26-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S26')

subplot(2,2,3)
rfplot(s_params_rebuilt,2,4); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,2,4) % submat of simulated S
hold off
% legend('off')
legend({'S24-rebuilt','S24-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S24')

subplot(2,2,4)
rfplot(s_params_rebuilt,2,8); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,2,8) % submat of simulated S
hold off
% legend('off')
legend({'S28-rebuilt','S28-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S28')

%% Line 34
figure('Name','Verify Rebuilt S params: Line 34')
subplot(2,2,1)
rfplot(s_params_rebuilt,3,3); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,3,3) % submat of simulated S
hold off
% legend('off')
legend({'S33-rebuilt','S33-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S33')

subplot(2,2,2)
rfplot(s_params_rebuilt,3,7); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,3,7) % submat of simulated S
hold off
% legend('off')
legend({'S37-rebuilt','S37-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S37')

subplot(2,2,3)
rfplot(s_params_rebuilt,3,4); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,3,4) % submat of simulated S
hold off
% legend('off')
legend({'S34-rebuilt','S34-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S34')

subplot(2,2,4)
rfplot(s_params_rebuilt,3,8); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params_1,3,8) % submat of simulated S
hold off
% legend('off')
legend({'S38-rebuilt','S38-simu'},'Location','best','NumColumns',1)
legend('boxoff')
title('S38')


