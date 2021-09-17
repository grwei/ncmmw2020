%% Basic information
% The verification of the following paper(Section II):
% M. K. Sampath, "On addressing the practical issues in the extraction of RLGC parameters for lossy multiconductor transmission lines using S-parameter models," 2008 IEEE-EPEP Electrical Performance of Electronic Packaging, San Jose, CA, 2008, pp. 259-262.
% Author Name: Wei Guorui
% Instructor Name: Xia Bin
% Created at: 2019-10-23 16:45
% Last modified: 2019-10-26 12:11
clc; clear; close all;
addpath(genpath('../code'))
addpath(genpath('../code/function'))

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
linelength = 0.015; % line length(meters).
filename_4line = 'data/four_line/d4Line_15mm_201910232329.s8p';
SE8P_Data = read(rfdata.data,filename_4line); % SingleEnded8PortData
freq = SE8P_Data.Freq;
freqpts = length(freq);
Z0 = SE8P_Data.Z0;
ZL = SE8P_Data.ZL;
ZS = SE8P_Data.ZS;
SE8P_Data.S_Parameters = snp2smp(SE8P_Data.S_Parameters,...
    Z0,[4 2 1 3 8 6 5 7]);
numLines = size(SE8P_Data.S_Parameters,1)/2;
s_params = sparameters(SE8P_Data);
S11_simu = squeeze(SE8P_Data.S_Parameters(1,1,:)); % used to find ouliers

%% Convert 2N¡Á2N S matrix to 2N¡Á2N T Matrix(ABCD) and Z matrix
T_mat = s2t_mod(SE8P_Data.S_Parameters, Z0); % notice the difference in "s2rlgc": S -> T.
Td = T_mat(numLines+1:end, numLines+1:end,:);
% Td(:,:,idx) = ((I-s11(:,:,idx))*(I+s22(:,:,idx))+s12(:,:,idx)*s21(:,:,idx))/(2*s21(:,:,idx));
Z_params = s2z(SE8P_Data.S_Parameters,Z0);

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
rlgc_params_matlab= s2rlgc_mod2(SE8P_Data.S_Parameters,linelength,freq,Z0);

%% fit
outliers = find(db(S11_simu) < -9.2);
[rlgc_params_fit] = rlgc_fit(rlgc_params,freq,outliers);

%% Rebuilt S params
% One-line method before fitted
s_params_rebuilt = snp2smp(sparameters(rlgc2s_mod(rlgc_params.R,rlgc_params.L,...
    rlgc_params.G,rlgc_params.C,linelength,freq,Z0),freq,Z0),[1 2 3 4 5 6 7 8],Z0);

% One-line method after fitted
s_params_fitted = snp2smp(sparameters(rlgc2s_mod(rlgc_params_fit.R,rlgc_params_fit.L,...
    rlgc_params_fit.G,rlgc_params_fit.C,linelength,freq,Z0),freq,Z0),[1 2 3 4 5 6 7 8],Z0);

% Compare with RF Toolbox("s2rlgc")
s_params_rebuilt_matlab = snp2smp(sparameters(rlgc2s_mod(rlgc_params_matlab.R,rlgc_params_matlab.L,...
    rlgc_params_matlab.G,rlgc_params_matlab.C,linelength,freq,Z0),freq,Z0),[1 2 3 4 5 6 7 8],Z0);

%% figure.

%% Simulated Scaterring matrix
figure('Name','Simulated [S] Matrix')
subplot(2,4,1)
rfplot(s_params,1,1);
hold on
rfplot(s_params,5,5);
rfplot(s_params,4,4);
rfplot(s_params,8,8);
hold off
legend({'S11','S55','S44','S88'},'Location','best','NumColumns',2)
legend('boxoff')
title('Reflection: S11,S55,S44,S88')

subplot(2,4,2)
rfplot(s_params,2,2);
hold on
rfplot(s_params,6,6);
rfplot(s_params,3,3);
rfplot(s_params,7,7);
hold off
legend({'S22','S66','S33','S77'},'Location','best','NumColumns',2)
legend('boxoff')
title('Reflection: S22,S66,S33,S77')

subplot(2,4,3)
rfplot(s_params,1,5);
hold on
rfplot(s_params,5,1);
rfplot(s_params,4,8);
rfplot(s_params,8,4);
hold off
legend({'S15','S51','S48','S84'},'Location','best','NumColumns',2)
legend('boxoff')
title('Insertion: S15,S51,S48,S84')

subplot(2,4,4)
rfplot(s_params,2,6);
hold on
rfplot(s_params,6,2);
rfplot(s_params,3,7);
rfplot(s_params,7,3);
hold off
legend({'S26','S62','S37','S73'},'Location','best','NumColumns',2)
legend('boxoff')
title('Insertion: S26,S62,S37,S73')

subplot(2,4,5)
rfplot(s_params,1,6);
hold on
rfplot(s_params,6,1);
rfplot(s_params,2,5);
rfplot(s_params,5,2);
hold off
legend({'S16','S61','S25','S52'},'Location','best','NumColumns',2)
legend('boxoff')
title('Coupling: S16,S61,S25,S52')

subplot(2,4,6)
rfplot(s_params,3,8);
hold on
rfplot(s_params,8,3);
rfplot(s_params,4,7);
rfplot(s_params,7,4);
hold off
legend({'S38','S83','S47','S74'},'Location','best','NumColumns',2)
legend('boxoff')
title('Coupling: S38,S83,S47,S74')

subplot(2,4,7)
rfplot(s_params,1,2);
hold on
rfplot(s_params,2,1);
rfplot(s_params,5,6);
rfplot(s_params,6,5);
hold off
legend({'S12','S21','S56','S65'},'Location','best','NumColumns',2)
legend('boxoff')
title('Isolation: S12,S21,S56,S65')

subplot(2,4,8)
rfplot(s_params,3,4);
hold on
rfplot(s_params,4,3);
rfplot(s_params,7,8);
rfplot(s_params,8,7);
hold off
legend({'S34','S43','S78','S87'},'Location','best','NumColumns',2)
legend('boxoff')
title('Isolation: S34,S43,S78,S87')

% figure. S params assumes to zero
figure('Name','Simulated [S] Matrix: Assumes to zero')
subplot(232)
rfplot(s_params,1,7);
hold on
rfplot(s_params,7,1);
rfplot(s_params,3,5);
rfplot(s_params,5,3);
rfplot(s_params,2,8);
rfplot(s_params,8,2);
rfplot(s_params,4,6);
rfplot(s_params,6,4);
hold off
legend({'S17','S71','S35','S53','S28','S82','S46','S64'},'Location','best','NumColumns',3)
legend('boxoff')
title("Mid Neighbor: S17,S71,S35,S53," + newline + "S28,S82,S46,S64")

subplot(235)
rfplot(s_params,1,3);
hold on
rfplot(s_params,3,1);
rfplot(s_params,2,4);
rfplot(s_params,4,2);
rfplot(s_params,5,7);
rfplot(s_params,7,5);
rfplot(s_params,6,8);
rfplot(s_params,8,6);
hold off
legend({'S13','S31','S24','S42','S57','S75','S68','S86'},'Location','best','NumColumns',3)
legend('boxoff')
title("Mid Neighbor: S13,S31,S24,S42," + newline + "S57,S75,S68,S86")

subplot(233)
rfplot(s_params,1,8);
hold on
rfplot(s_params,8,1);
rfplot(s_params,4,5);
rfplot(s_params,5,4);
hold off
legend({'S18','S81','S45','S54'},'Location','best','NumColumns',2)
legend('boxoff')
title("Far Neighbor: S18,S81,S45,S54")

subplot(236)
rfplot(s_params,1,4);
hold on
rfplot(s_params,4,1);
rfplot(s_params,5,8);
rfplot(s_params,8,5);
hold off
legend({'S14','S41','S58','S85'},'Location','best','NumColumns',2)
legend('boxoff')
title("Far Neighbor: S14,S41,S58,S85")

subplot(231)
rfplot(s_params,2,7);
hold on
rfplot(s_params,7,2);
rfplot(s_params,3,6);
rfplot(s_params,6,3);
hold off
legend({'S27','S72','S36','S63'},'Location','best','NumColumns',2)
legend('boxoff')
title('Close Neighbor: S27,S72,S36,S63')

subplot(234)
rfplot(s_params,2,3);
hold on
rfplot(s_params,3,2);
rfplot(s_params,6,7);
rfplot(s_params,7,6);
hold off
legend({'S23','S32','S67','S76'},'Location','best','NumColumns',2)
legend('boxoff')
title('Close Neighbor: S23,S32,S67,S76')

%% Extracted RLGC params: One-line method
figure('Name','R Matrix: One-line method')
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

figure('Name','L Matrix: One-line method')
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

figure('Name','C Matrix: One-line method')
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

figure('Name','G Matrix: One-line method')
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

%% 
%% Verify Rebuilt S params(One line method)

%% Line 12
figure('Name','Verify Rebuilt S params: Line 12')
subplot(2,2,1)
rfplot(s_params_rebuilt,1,1); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,1) % submat of simulated S
rfplot(s_params_fitted,1,1); 
hold off
% legend('off')
legend({'S11-rebuilt','S11-simu','S11-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S11')

subplot(2,2,2)
rfplot(s_params_rebuilt,1,2); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,2) % submat of simulated S
rfplot(s_params_fitted,1,2); 
hold off
% legend('off')
legend({'S12-rebuilt','S12-simu','S12-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S12')

subplot(2,2,3)
rfplot(s_params_rebuilt,1,5); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,5) % submat of simulated S
rfplot(s_params_fitted,1,5); 
hold off
% legend('off')
legend({'S15-rebuilt','S15-simu','S15-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S15')

subplot(2,2,4)
rfplot(s_params_rebuilt,1,6); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,6) % submat of simulated S
rfplot(s_params_fitted,1,6); 
hold off
% legend('off')
legend({'S16-rebuilt','S16-simu','S16-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S16')

%% Line 13
figure('Name','Verify Rebuilt S params: Line 13')
subplot(2,2,1)
rfplot(s_params_rebuilt,1,1); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,1) % submat of simulated S
rfplot(s_params_fitted,1,1); 
hold off
% legend('off')
legend({'S11-rebuilt','S11-simu','S11-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S11')

subplot(2,2,2)
rfplot(s_params_rebuilt,1,3); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,3) % submat of simulated S
rfplot(s_params_fitted,1,3); 
hold off
% legend('off')
legend({'S13-rebuilt','S13-simu','S13-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S13')

subplot(2,2,3)
rfplot(s_params_rebuilt,1,5); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,5) % submat of simulated S
rfplot(s_params_fitted,1,5); 
hold off
% legend('off')
legend({'S15-rebuilt','S15-simu','S15-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S15')

subplot(2,2,4)
rfplot(s_params_rebuilt,1,7); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,7) % submat of simulated S
rfplot(s_params_fitted,1,7); 
hold off
% legend('off')
legend({'S17-rebuilt','S17-simu','S17-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S17')

%% Line 14
figure('Name','Verify Rebuilt S params: Line 14')
subplot(2,2,1)
rfplot(s_params_rebuilt,1,1); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,1) % submat of simulated S
rfplot(s_params_fitted,1,1); 
hold off
% legend('off')
legend({'S11-rebuilt','S11-simu','S11-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S11')

subplot(2,2,2)
rfplot(s_params_rebuilt,1,4); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,4) % submat of simulated S
rfplot(s_params_fitted,1,4); 
hold off
% legend('off')
legend({'S14-rebuilt','S14-simu','S14-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S14')

subplot(2,2,3)
rfplot(s_params_rebuilt,1,5); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,5) % submat of simulated S
rfplot(s_params_fitted,1,5); 
hold off
% legend('off')
legend({'S15-rebuilt','S15-simu','S15-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S15')

subplot(2,2,4)
rfplot(s_params_rebuilt,1,8); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,1,8) % submat of simulated S
rfplot(s_params_fitted,1,8); 
hold off
% legend('off')
legend({'S18-rebuilt','S18-simu','S18-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S18')

%% Line 23
figure('Name','Verify Rebuilt S params: Line 23')
subplot(2,2,1)
rfplot(s_params_rebuilt,2,2); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,2,2) % submat of simulated S
rfplot(s_params_fitted,2,2); 
hold off
% legend('off')
legend({'S22-rebuilt','S22-simu','S22-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S22')

subplot(2,2,2)
rfplot(s_params_rebuilt,2,6); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,2,6) % submat of simulated S
rfplot(s_params_fitted,2,6); 
hold off
% legend('off')
legend({'S26-rebuilt','S26-simu','S26-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S26')

subplot(2,2,3)
rfplot(s_params_rebuilt,2,3); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,2,3) % submat of simulated S
rfplot(s_params_fitted,2,3); 
hold off
% legend('off')
legend({'S23-rebuilt','S23-simu','S23-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S23')

subplot(2,2,4)
rfplot(s_params_rebuilt,2,7); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,2,7) % submat of simulated S
rfplot(s_params_fitted,2,7); 
hold off
% legend('off')
legend({'S27-rebuilt','S27-simu','S27-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S27')

%% Line 24
figure('Name','Verify Rebuilt S params: Line 24')
subplot(2,2,1)
rfplot(s_params_rebuilt,2,2); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,2,2) % submat of simulated S
rfplot(s_params_fitted,2,2); 
hold off
% legend('off')
legend({'S22-rebuilt','S22-simu','S22-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S22')

subplot(2,2,2)
rfplot(s_params_rebuilt,2,6); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,2,6) % submat of simulated S
rfplot(s_params_fitted,2,6); 
hold off
% legend('off')
legend({'S26-rebuilt','S26-simu','S26-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S26')

subplot(2,2,3)
rfplot(s_params_rebuilt,2,4); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,2,4) % submat of simulated S
rfplot(s_params_fitted,2,4); 
hold off
% legend('off')
legend({'S24-rebuilt','S24-simu','S24-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S24')

subplot(2,2,4)
rfplot(s_params_rebuilt,2,8); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,2,8) % submat of simulated S
rfplot(s_params_fitted,2,8); 
hold off
% legend('off')
legend({'S28-rebuilt','S28-simu','S28-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S28')

%% Line 34
figure('Name','Verify Rebuilt S params: Line 34')
subplot(2,2,1)
rfplot(s_params_rebuilt,3,3); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,3,3) % submat of simulated S
rfplot(s_params_fitted,3,3); 
hold off
% legend('off')
legend({'S33-rebuilt','S33-simu','S33-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S33')

subplot(2,2,2)
rfplot(s_params_rebuilt,3,7); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,3,7) % submat of simulated S
rfplot(s_params_fitted,3,7); 
hold off
% legend('off')
legend({'S37-rebuilt','S37-simu','S37-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S37')

subplot(2,2,3)
rfplot(s_params_rebuilt,3,4); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,3,4) % submat of simulated S
rfplot(s_params_fitted,3,4); 
hold off
% legend('off')
legend({'S34-rebuilt','S34-simu','S34-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S34')

subplot(2,2,4)
rfplot(s_params_rebuilt,3,8); % submat of rebuilt S(from total rlgc)
hold on
rfplot(s_params,3,8) % submat of simulated S
rfplot(s_params_fitted,3,8); 
hold off
% legend('off')
legend({'S38-rebuilt','S38-simu','S38-fitted'},'Location','best','NumColumns',1)
legend('boxoff')
title('S38')

%% Extracted RLGC params: One-line method(fitted)
figure('Name','R Matrix: One-line method(fitted)')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params_fit.R(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.R(2,2,:)))
plot(freq,squeeze(rlgc_params_fit.R(3,3,:)))
plot(freq,squeeze(rlgc_params_fit.R(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R11,R22,R33,R44');
legend({'R11','R22','R33','R44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params_fit.R(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.R(2,1,:)))
plot(freq,squeeze(rlgc_params_fit.R(3,4,:)))
plot(freq,squeeze(rlgc_params_fit.R(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R12,R21,R34,R43');
legend({'R12','R21','R34','R43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params_fit.R(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.R(3,1,:)))
plot(freq,squeeze(rlgc_params_fit.R(2,4,:)))
plot(freq,squeeze(rlgc_params_fit.R(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R13,R31,R24,R42');
legend({'R13','R31','R24','R42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params_fit.R(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.R(4,1,:)))
plot(freq,squeeze(rlgc_params_fit.R(2,3,:)))
plot(freq,squeeze(rlgc_params_fit.R(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R_{ij}(Ohm/m)');
title('R14,R41,R23,R32');
legend({'R14','R41','R23','R32'},'Location','best','NumColumns',2)
legend('boxoff')

figure('Name','L Matrix: One-line method(fitted)')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params_fit.L(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.L(2,2,:)))
plot(freq,squeeze(rlgc_params_fit.L(3,3,:)))
plot(freq,squeeze(rlgc_params_fit.L(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L11,L22,L33,L44');
legend({'L11','L22','L33','L44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params_fit.L(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.L(2,1,:)))
plot(freq,squeeze(rlgc_params_fit.L(3,4,:)))
plot(freq,squeeze(rlgc_params_fit.L(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L12,L21,L34,L43');
legend({'L12','L21','L34','L43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params_fit.L(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.L(3,1,:)))
plot(freq,squeeze(rlgc_params_fit.L(2,4,:)))
plot(freq,squeeze(rlgc_params_fit.L(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L13,L31,L24,L42');
legend({'L13','L31','L24','L42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params_fit.L(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.L(4,1,:)))
plot(freq,squeeze(rlgc_params_fit.L(2,3,:)))
plot(freq,squeeze(rlgc_params_fit.L(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L_{ij}(H/m)');
title('L14,L41,L23,L32');
legend({'L14','L41','L23','L32'},'Location','best','NumColumns',2)
legend('boxoff')

figure('Name','C Matrix: One-line method(fitted)')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params_fit.C(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.C(2,2,:)))
plot(freq,squeeze(rlgc_params_fit.C(3,3,:)))
plot(freq,squeeze(rlgc_params_fit.C(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C11,C22,C33,C44');
legend({'C11','C22','C33','C44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params_fit.C(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.C(2,1,:)))
plot(freq,squeeze(rlgc_params_fit.C(3,4,:)))
plot(freq,squeeze(rlgc_params_fit.C(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C12,C21,C34,C43');
legend({'C12','C21','C34','C43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params_fit.C(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.C(3,1,:)))
plot(freq,squeeze(rlgc_params_fit.C(2,4,:)))
plot(freq,squeeze(rlgc_params_fit.C(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C13,C31,C24,C42');
legend({'C13','C31','C24','C42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params_fit.C(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.C(4,1,:)))
plot(freq,squeeze(rlgc_params_fit.C(2,3,:)))
plot(freq,squeeze(rlgc_params_fit.C(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C_{ij}(F/m)');
title('C14,C41,C23,C32');
legend({'C14','C41','C23','C32'},'Location','best','NumColumns',2)
legend('boxoff')

figure('Name','G Matrix: One-line method(fitted)')
subplot(2,2,1)
plot(freq,squeeze(rlgc_params_fit.G(1,1,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.G(2,2,:)))
plot(freq,squeeze(rlgc_params_fit.G(3,3,:)))
plot(freq,squeeze(rlgc_params_fit.G(4,4,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G11,G22,G33,G44');
legend({'G11','G22','G33','G44'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,2)
plot(freq,squeeze(rlgc_params_fit.G(1,2,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.G(2,1,:)))
plot(freq,squeeze(rlgc_params_fit.G(3,4,:)))
plot(freq,squeeze(rlgc_params_fit.G(4,3,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G12,G21,G34,G43');
legend({'G12','G21','G34','G43'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,3)
plot(freq,squeeze(rlgc_params_fit.G(1,3,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.G(3,1,:)))
plot(freq,squeeze(rlgc_params_fit.G(2,4,:)))
plot(freq,squeeze(rlgc_params_fit.G(4,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G13,G31,G24,G42');
legend({'G13','G31','G24','G42'},'Location','best','NumColumns',2)
legend('boxoff')

subplot(2,2,4)
plot(freq,squeeze(rlgc_params_fit.G(1,4,:)))
hold on
plot(freq,squeeze(rlgc_params_fit.G(4,1,:)))
plot(freq,squeeze(rlgc_params_fit.G(2,3,:)))
plot(freq,squeeze(rlgc_params_fit.G(3,2,:)))
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G_{ij}(S/m)');
title('G14,G41,G23,G32');
legend({'G14','G41','G23','G32'},'Location','best','NumColumns',2)
legend('boxoff')

%% Verify fitted S params
%% Simulated Scaterring matrix
figure('Name','Fitted [S] Matrix')
subplot(2,4,1)
rfplot(s_params,1,1);
hold on
rfplot(s_params,5,5);
rfplot(s_params,4,4);
rfplot(s_params,8,8);
rfplot(s_params_fitted,1,1);
rfplot(s_params_fitted,5,5);
rfplot(s_params_fitted,4,4);
rfplot(s_params_fitted,8,8);
hold off
legend({'S11','S55','S44','S88','S11-fit','S55-fit','S44-fit','S88-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title('Reflection: S11,S55,S44,S88')

subplot(2,4,2)
rfplot(s_params,2,2);
hold on
rfplot(s_params,6,6);
rfplot(s_params,3,3);
rfplot(s_params,7,7);
rfplot(s_params_fitted,2,2);
rfplot(s_params_fitted,6,6);
rfplot(s_params_fitted,3,3);
rfplot(s_params_fitted,7,7);
hold off
legend({'S22','S66','S33','S77','S22-fit','S66-fit','S33-fit','S77-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title('Reflection: S22,S66,S33,S77')

subplot(2,4,3)
rfplot(s_params,1,5);
hold on
rfplot(s_params,5,1);
rfplot(s_params,4,8);
rfplot(s_params,8,4);
rfplot(s_params_fitted,1,5);
rfplot(s_params_fitted,5,1);
rfplot(s_params_fitted,4,8);
rfplot(s_params_fitted,8,4);
hold off
legend({'S15','S51','S48','S84','S15-fit','S51-fit','S48-fit','S84-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title('Insertion: S15,S51,S48,S84')

subplot(2,4,4)
rfplot(s_params,2,6);
hold on
rfplot(s_params,6,2);
rfplot(s_params,3,7);
rfplot(s_params,7,3);
rfplot(s_params_fitted,2,6);
rfplot(s_params_fitted,6,2);
rfplot(s_params_fitted,3,7);
rfplot(s_params_fitted,7,3);
hold off
legend({'S26','S62','S37','S73','S26-fit','S62-fit','S37-fit','S73-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title('Insertion: S26,S62,S37,S73')

subplot(2,4,5)
rfplot(s_params,1,6);
hold on
rfplot(s_params,6,1);
rfplot(s_params,2,5);
rfplot(s_params,5,2);
rfplot(s_params_fitted,1,6);
rfplot(s_params_fitted,6,1);
rfplot(s_params_fitted,2,5);
rfplot(s_params_fitted,5,2);
hold off
legend({'S16','S61','S25','S52','S16-fit','S61-fit','S25-fit','S52-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title('Coupling: S16,S61,S25,S52')

subplot(2,4,6)
rfplot(s_params,3,8);
hold on
rfplot(s_params,8,3);
rfplot(s_params,4,7);
rfplot(s_params,7,4);
rfplot(s_params_fitted,3,8);
rfplot(s_params_fitted,8,3);
rfplot(s_params_fitted,4,7);
rfplot(s_params_fitted,7,4);
hold off
legend({'S38','S83','S47','S74','S38-fit','S83-fit','S47-fit','S74-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title('Coupling: S38,S83,S47,S74')

subplot(2,4,7)
rfplot(s_params,1,2);
hold on
rfplot(s_params,2,1);
rfplot(s_params,5,6);
rfplot(s_params,6,5);
rfplot(s_params_fitted,1,2);
rfplot(s_params_fitted,2,1);
rfplot(s_params_fitted,5,6);
rfplot(s_params_fitted,6,5);
hold off
legend({'S12','S21','S56','S65','S12-fit','S21-fit','S56-fit','S65-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title('Isolation: S12,S21,S56,S65')

subplot(2,4,8)
rfplot(s_params,3,4);
hold on
rfplot(s_params,4,3);
rfplot(s_params,7,8);
rfplot(s_params,8,7);
rfplot(s_params_fitted,3,4);
rfplot(s_params_fitted,4,3);
rfplot(s_params_fitted,7,8);
rfplot(s_params_fitted,8,7);
hold off
legend({'S34','S43','S78','S87','S34-fit','S43-fit','S78-fit','S87-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title('Isolation: S34,S43,S78,S87')

% figure. S params assumes to zero
figure('Name','Fitted [S] Matrix: Assumes to zero')
subplot(232)
rfplot(s_params,1,7);
hold on
rfplot(s_params,7,1);
rfplot(s_params,3,5);
rfplot(s_params,5,3);
rfplot(s_params,2,8);
rfplot(s_params,8,2);
rfplot(s_params,4,6);
rfplot(s_params,6,4);
rfplot(s_params_fitted,1,7);
rfplot(s_params_fitted,7,1);
rfplot(s_params_fitted,3,5);
rfplot(s_params_fitted,5,3);
rfplot(s_params_fitted,2,8);
rfplot(s_params_fitted,8,2);
rfplot(s_params_fitted,4,6);
rfplot(s_params_fitted,6,4);
hold off
legend({'S17','S71','S35','S53','S28','S82','S46','S64','S17-fit','S71-fit','S35-fit','S53-fit','S28-fit','S82-fit','S46-fit','S64-fit'},'Location','best','NumColumns',4)
legend('boxoff')
title("Mid Neighbor: S17,S71,S35,S53," + newline + "S28,S82,S46,S64")

subplot(235)
rfplot(s_params,1,3);
hold on
rfplot(s_params,3,1);
rfplot(s_params,2,4);
rfplot(s_params,4,2);
rfplot(s_params,5,7);
rfplot(s_params,7,5);
rfplot(s_params,6,8);
rfplot(s_params,8,6);
rfplot(s_params_fitted,1,3);
rfplot(s_params_fitted,3,1);
rfplot(s_params_fitted,2,4);
rfplot(s_params_fitted,4,2);
rfplot(s_params_fitted,5,7);
rfplot(s_params_fitted,7,5);
rfplot(s_params_fitted,6,8);
rfplot(s_params_fitted,8,6);
hold off
legend({'S13','S31','S24','S42','S57','S75','S68','S86','S13-fit','S31-fit','S24-fit','S42-fit','S57-fit','S75-fit','S68-fit','S86-fit'},'Location','best','NumColumns',4)
legend('boxoff')
title("Mid Neighbor: S13,S31,S24,S42," + newline + "S57,S75,S68,S86")

subplot(233)
rfplot(s_params,1,8);
hold on
rfplot(s_params,8,1);
rfplot(s_params,4,5);
rfplot(s_params,5,4);
rfplot(s_params_fitted,1,8);
rfplot(s_params_fitted,8,1);
rfplot(s_params_fitted,4,5);
rfplot(s_params_fitted,5,4);
hold off
legend({'S18','S81','S45','S54','S18-fit','S81-fit','S45-fit','S54-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title("Far Neighbor: S18,S81,S45,S54")

subplot(236)
rfplot(s_params,1,4);
hold on
rfplot(s_params,4,1);
rfplot(s_params,5,8);
rfplot(s_params,8,5);
rfplot(s_params_fitted,1,4);
rfplot(s_params_fitted,4,1);
rfplot(s_params_fitted,5,8);
rfplot(s_params_fitted,8,5);
hold off
legend({'S14','S41','S58','S85','S14-fit','S41-fit','S58-fit','S85-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title("Far Neighbor: S14,S41,S58,S85")

subplot(231)
rfplot(s_params,2,7);
hold on
rfplot(s_params,7,2);
rfplot(s_params,3,6);
rfplot(s_params,6,3);
rfplot(s_params_fitted,2,7);
rfplot(s_params_fitted,7,2);
rfplot(s_params_fitted,3,6);
rfplot(s_params_fitted,6,3);
hold off
legend({'S27','S72','S36','S63','S27-fit','S72-fit','S36-fit','S63-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title('Close Neighbor: S27,S72,S36,S63')

subplot(234)
rfplot(s_params,2,3);
hold on
rfplot(s_params,3,2);
rfplot(s_params,6,7);
rfplot(s_params,7,6);
rfplot(s_params_fitted,2,3);
rfplot(s_params_fitted,3,2);
rfplot(s_params_fitted,6,7);
rfplot(s_params_fitted,7,6);
hold off
legend({'S23','S32','S67','S76','S23-fit','S32-fit','S67-fit','S76-fit'},'Location','best','NumColumns',2)
legend('boxoff')
title('Close Neighbor: S23,S32,S67,S76')