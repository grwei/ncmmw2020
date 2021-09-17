%% Basic information
% Extract RLGC of coupled TL from S using Mixed-Mode ABCD method
% Author Name: Wei Guorui
% Instructor Name: Xia Bin
% Created in: 2020-01-01 11:16
% Last modified: 2020-01-01 hh:mm
clc; clear; close all;
addpath(genpath('../code'))

%% Port numbering
% Single-Ended
% 	           PORT NUMBERING (Classic style)
% 	    PORT#                               PORT#
% 	    Z1 = Z3 = 50.00                     Z2 = Z4 = 50.00
% 	           --------------------------
% 	          |                          |
% 	      1 --|                          |-- 2
% 	          |                          |
% 	      3 --|                          |-- 4
% 	          |                          |
% 	           --------------------------

% Mixed Mode S-Parameters
% 	    Zdiff = 100.00                      Zcomm = 25.00
% 	    PORT#								PORT#
% 	           --------------------------
% 	          |                          |
% 	        --|                          |--
% 	     1    |                          |     2
% 	        --|                          |--
% 	          |                          |
% 	           --------------------------

%% Initialize
% Load Scattering matrix
linelength = 0.018; % line length(meters).
filename_2line = 'data/SI/d20200101.s4p';
SE4P_Data = read(rfdata.data,filename_2line); % SingleEnded4PortData
freq = SE4P_Data.Freq;
freqpts = length(freq);
Z0 = SE4P_Data.Z0; % Reference impedance
ZL = SE4P_Data.ZL;
ZS = SE4P_Data.ZS;
SE4P_Data.S_Parameters = snp2smp(SE4P_Data.S_Parameters,...
    Z0,[1 2 3 4]); % Single-Ended Scattering parameters
numLines = size(SE4P_Data.S_Parameters,1)/2; % (nunLines + 1) Conductors
s_params = sparameters(SE4P_Data);

%% Extract RLCG using Mixed-Mode ABCD method
% 1. Convert Single-Ended S to Mixed-Mode S
[smm.dd,smm.dc,smm.cd,smm.cc] = ...
    s2smm(snp2smp(s_params.Parameters,Z0,[1 2 3 4],Z0),2);

% 2. Extract \gamma and Zc for each mode
% 3. Combined odd/even mode to get total RLCG
[rlgc_params] = mms2rlgc(smm,linelength,freq,Z0);

%% Validate
% Load RLCG extracted by Cadence Sigrity PowerSI
% allocate memory
rlgc_PowerSI.R = zeros(size(rlgc_params.R));
rlgc_PowerSI.L = zeros(size(rlgc_params.L));
rlgc_PowerSI.C = zeros(size(rlgc_params.C));
rlgc_PowerSI.G = zeros(size(rlgc_params.G));
% load data
filename_PowerSI = 'data/SI/Transmission_RLGC_res_20200101.csv';
opts = detectImportOptions(filename_PowerSI);
rlgc_PowerSI_mat = readtable(filename_PowerSI);
for freqIdx = 1:freqpts
    rlgc_PowerSI.R(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,3}/linelength;
    rlgc_PowerSI.R(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,4}/linelength;
    rlgc_PowerSI.R(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,4}/linelength;
    rlgc_PowerSI.R(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-3,5}/linelength;
    rlgc_PowerSI.L(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,3}/linelength;
    rlgc_PowerSI.L(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,4}/linelength;
    rlgc_PowerSI.L(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,4}/linelength;
    rlgc_PowerSI.L(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-2,5}/linelength;    
    rlgc_PowerSI.G(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,3}/linelength;
    rlgc_PowerSI.G(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,4}/linelength;
    rlgc_PowerSI.G(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,4}/linelength;
    rlgc_PowerSI.G(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-1,5}/linelength;
    rlgc_PowerSI.C(1,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,3}/linelength;
    rlgc_PowerSI.C(1,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,4}/linelength;
    rlgc_PowerSI.C(2,1,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,4}/linelength;
    rlgc_PowerSI.C(2,2,freqIdx) = rlgc_PowerSI_mat{4*freqIdx-0,5}/linelength;
end

%% Figure
%% 1. Simulated S parameters
figure('Name','Simulated [S] Matrix')
sgtitle('Simulated S parameters using Polar Si9000')
subplot(2,2,1)
rfplot(s_params,1,1);
hold on
rfplot(s_params,2,2);
rfplot(s_params,3,3);
rfplot(s_params,4,4);
hold off
legend({'S11','S22','S33','S44'},'Location','best','NumColumns',2)
legend('boxoff')
title('S11,S22,S33,S44')

subplot(2,2,2)
rfplot(s_params,1,2);
hold on
rfplot(s_params,2,1);
rfplot(s_params,3,4);
rfplot(s_params,4,3);
hold off
legend({'S12','S21','S34','S43'},'Location','best','NumColumns',2)
legend('boxoff')
title('S12,S21,S34,S43')

subplot(2,2,3)
rfplot(s_params,1,3);
hold on
rfplot(s_params,3,1);
rfplot(s_params,2,4);
rfplot(s_params,4,2);
hold off
legend({'S13','S31','S24','S42'},'Location','best','NumColumns',2)
legend('boxoff')
title('S13,S31,S24,S42')

subplot(2,2,4)
rfplot(s_params,1,4);
hold on
rfplot(s_params,4,1);
rfplot(s_params,2,3);
rfplot(s_params,3,2);
hold off
legend({'S14','S41','S23','S32'},'Location','best','NumColumns',2)
legend('boxoff')
title('S14,S41,S23,S32')

%% 2. Propagation constants and characteristic impedances of even/odd mode
figure('Name','Propagation constants and characteristic impedances of even/odd mode')
sgtitle({'Propagation constants and characteristic impedances','of even/odd mode'})

subplot(2,2,1)
plot(freq/1e9,real(rlgc_params.even.alpha),'-')
hold on
plot(freq/1e9,real(rlgc_params.odd.alpha),'--')
grid on
xlabel('Freq(GHz)');
ylabel('\alpha(Np/m)');
title('\alpha:even/odd mode');
legend({'\alpha_e','\alpha_o'},'Location','best')

subplot(2,2,2)
plot(freq/1e9,real(rlgc_params.even.beta),'k-')
hold on
plot(freq/1e9,real(rlgc_params.odd.beta),'g--')
grid on
xlabel('Freq(GHz)');
ylabel('\beta(rad/m)');
title('\beta:even/odd mode');
legend({'\beta_e','\beta_o'},'Location','best')

subplot(2,2,3)
plot(freq/1e9,real(rlgc_params.even.Zc))
hold on
plot(freq/1e9,real(rlgc_params.odd.Zc))
grid on
xlabel('Freq(GHz)');
ylabel('\Re(Zc)/Ohms');
title('\Re(Zc):even/odd mode');
legend({'\Re(Zc,e)','\Re(Zc,o)'},'Location','best')

subplot(2,2,4)
plot(freq/1e9,imag(rlgc_params.even.Zc))
hold on
plot(freq/1e9,imag(rlgc_params.odd.Zc))
grid on
xlabel('Freq(GHz)');
ylabel('\Im(Zc)/Ohms');
title('\Im(Zc):even/odd mode');
legend({'\Im(Zc,e)','\Im(Zc,o)'},'Location','best')

%% 3. RLGC compared with Cadence PowerSI
% R, L
figure('Name','R, L matrix')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: R, L Matrix'})
subplot(221)
plot(freq/1e9,squeeze(rlgc_params.R(1,1,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.R(1,1,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('R11(Ohms/m)');
title('R11 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(222)
plot(freq/1e9,squeeze(rlgc_params.R(1,2,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.R(1,2,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('R12(Ohms/m)');
title('R12 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(223)
plot(freq/1e9,squeeze(rlgc_params.L(1,1,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.L(1,1,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('L11(H/m)');
title('L11 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(224)
plot(freq/1e9,squeeze(rlgc_params.L(1,2,:)),'k-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.L(1,2,:)),'g--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('L12(H/m)');
title('L12 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

% C, G
figure('Name','C, G matrix')
sgtitle({'Comparison Between Proposed Algorithm and';' Cadence Sigrity PowerSI: C, G Matrix'})
subplot(221)
plot(freq/1e9,squeeze(rlgc_params.C(1,1,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.C(1,1,:)),'k--')
hold off
grid on
ave_y = mean(squeeze(rlgc_params.C(1,1,:)));
del_r = 4e-5;
ylim([ave_y*(1-del_r), ave_y*(1+del_r)])
xlabel('Freq(GHz)');
ylabel('C11(F/m)');
title('C11 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(222)
plot(freq/1e9,squeeze(rlgc_params.C(1,2,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.C(1,2,:)),'k--')
hold off
grid on
ave_y = mean(squeeze(rlgc_params.C(1,2,:)));
del_r = 5e-4;
ylim([ave_y*(1+del_r), ave_y*(1-del_r)])
xlabel('Freq(GHz)');
ylabel('C12(F/m)');
title('C12 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(223)
plot(freq/1e9,squeeze(rlgc_params.G(1,1,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.G(1,1,:)),'k--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('G11(S/m)');
title('G11 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')

subplot(224)
plot(freq/1e9,squeeze(rlgc_params.G(1,2,:)),'m-')
hold on
plot(freq/1e9,squeeze(rlgc_PowerSI.G(1,2,:)),'k--')
hold off
grid on
xlabel('Freq(GHz)');
ylabel('G12(S/m)');
title('G12 Comparison');
legend({'Proposed Algorithm','Cadence Sigrity PowerSI'},'Location','best','NumColumns',1)
legend('boxoff')
