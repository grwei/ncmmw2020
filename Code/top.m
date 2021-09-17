%% Basic information
% Given the 4-Port S parameters and port characteristic impedance of a CPL network,
% extract the spice-Parameters:tabular RLGC model(W-element)
% Author Name: Wei Guorui
% Instructor Name: Xia Bin
% Created in: 2019-07-16 22:53
% Last modified: 2019-07-18 22:17
clc; clear; close all;
addpath(genpath('../code'))

%% Initialize
LL = 0.0254; % the length of the transmission line in meters.
filename_SI = '/data/SI/d07161350.xlsb';
[status,sheets_SI,xlFormat] = xlsfinfo(filename_SI);
data_SI_MMS = readtable(filename_SI,'sheet',string(sheets_SI(1))); % Mixed Mode S-Parameters
data_SI_4PS = readtable(filename_SI,'sheet',string(sheets_SI(2))); % 4 Port S-Parameters
data_SI_OM = readtable(filename_SI,'sheet',string(sheets_SI(3))); % Odd Mode
data_SI_EM = readtable(filename_SI,'sheet',string(sheets_SI(4))); % Even Mode
data_SI_SP = readtable(filename_SI,'sheet',string(sheets_SI(5))); % spice-Parameters:tabular RLGC model(W-element)
% opts = detectImportOptions(filename_SI,'sheet',string(sheets_SI(5)));
% freq =  data_SI_SP.Frequency_Hz_;

%% Load the simulated data(by Si9000)
freq =  data_SI_SP{:,1}; % Specify the vector of M frequencies over which the S-parameter array s_params is defined.
Z0 = 50; % Ports' reference(characteristic?) impedance in ohms
Zdiff = Z0 * 2; % ?
Zcomm = Z0 / 2; % ?
%% 4 Port S-Parameters
% 	           PORT NUMBERING (Modern style)
% 	    PORT#                               PORT#
% 	    Z1 = Z3 = 50.00                     Z2 = Z4 = 50.00
% 	           --------------------------
% 	          |                          |
% 	      1 --|                          |-- 2
% 	          |                          |
% 	      3 --|                          |-- 4
% 	          |                          |
% 	           --------------------------
S11 = complex(data_SI_4PS{:,2},data_SI_4PS{:,3});
S12 = complex(data_SI_4PS{:,4},data_SI_4PS{:,5});
S13 = complex(data_SI_4PS{:,6},data_SI_4PS{:,7});
S14 = complex(data_SI_4PS{:,8},data_SI_4PS{:,9});
% Scattering parameters of transmission line: 2N-by-2N-by-M array
s_params(4,4,:) = S11;s_params(3,3,:) = S11;s_params(2,2,:) = S11;s_params(1,1,:) = S11;
s_params(1,2,:) = S12;s_params(2,1,:) = S12;s_params(3,4,:) = S12;s_params(4,3,:) = S12;
s_params(1,3,:) = S13;s_params(3,1,:) = S13;s_params(2,4,:) = S13;s_params(4,2,:) = S13;
s_params(1,4,:) = S14;s_params(4,1,:) = S14;s_params(2,3,:) = S14;s_params(3,2,:) = S14;

%% Mixed Mode S-Parameters
% 	    Zdiff = 100.00                      Zcomm = 25.00
% 	    PORT#								PORT#
% 	           --------------------------
% 	          |                          |
% 	        --|                          |--   
% 	     1    |                          |     2
% 	        --|                          |--  
% 	          |                          |
% 	           --------------------------
SDD11 = complex(data_SI_MMS{:,2},data_SI_MMS{:,3});
SDD12 = complex(data_SI_MMS{:,4},data_SI_MMS{:,5});
SCC11 = complex(data_SI_MMS{:,26},data_SI_MMS{:,27});
SCC12 = complex(data_SI_MMS{:,28},data_SI_MMS{:,29});


%% spice-Parameters:tabular RLGC model(W-element)
R11_sim = data_SI_SP{:,2}; % tabular RLGC model(W-element)
R12_sim = data_SI_SP{:,3};
L11_sim = data_SI_SP{:,6};
L12_sim = data_SI_SP{:,7};
G11_sim = data_SI_SP{:,10};
G12_sim = -data_SI_SP{:,11};
C11_sim = data_SI_SP{:,14};
C12_sim = -data_SI_SP{:,15};

%% Coupled-Transmission Line Models
% 0. Test: Using RF Toolbox function "s2rlgc" & "rlgc2s"
rlgc_params_matlab = my_s2rlgc(snp2smp(s_params,Z0,[1 3 2 4]),LL,freq,Z0); % test toolbox function "s2rlgc"
s_params_matlab = snp2smp(rlgc2s(rlgc_params_matlab.R,rlgc_params_matlab.L,rlgc_params_matlab.G,rlgc_params_matlab.C,LL,freq,Z0),Z0,[1 3 2 4]);
% 1. Convert single-ended S-parameters to mixed-mode S-parameters
s_cc(2,2,:) = s_params(1,1,:) + s_params(1,3,:);s_cc(1,1,:) = s_cc(2,2,:);
s_cc(1,2,:) = s_params(1,2,:) + s_params(1,4,:);s_cc(2,1,:) = s_cc(1,2,:);
s_dd(2,2,:) = s_params(1,1,:) - s_params(1,3,:);s_dd(1,1,:) = s_dd(2,2,:);
s_dd(1,2,:) = s_params(1,2,:) - s_params(1,4,:);s_dd(2,1,:) = s_dd(1,2,:);
% what about s_cd? =0.
s_cd = zeros(size(s_cc));
s_dc = zeros(size(s_dd));
% verify using function"s2smm"
[vr_s_dd,vr_s_dc,vr_s_cd,vr_s_cc] = s2smm(s_params,1); % Convert single-ended S-parameters to mixed-mode S-parameters
chk_s_dd11 = min(vr_s_dd(1,1,:) == s_dd(1,1,:)); % chk == 1 means all are equal
chk_s_cc11 = min(vr_s_cc(1,1,:) == s_cc(1,1,:)); % chk == 1 means all are equal
chk_s_dd12 = min(vr_s_dd(1,2,:) == s_dd(1,2,:)); % chk == 1 means all are equal
chk_s_cc12 = min(vr_s_cc(1,2,:) == s_cc(1,2,:)); % chk == 1 means all are equal
chk_s_cd11 = min(vr_s_cd(1,1,:) == s_cd(1,1,:)); % chk == 1 means all are equal
chk_s_dc11 = min(vr_s_dc(1,1,:) == s_dc(1,1,:)); % chk == 1 means all are equal
chk_s_cd12 = min(vr_s_cd(1,2,:) == s_cd(1,2,:)); % chk == 1 means all are equal
chk_s_dc12 = min(vr_s_dc(1,2,:) == s_dc(1,2,:)); % chk == 1 means all are equal

% 2. Extracting the RLCG parameters for both the Diff-Diff mode and the
% Com-Com mode.
% rlgc_params = s2rlgc(s_cc,LL,freq,Zcomm); % 官方RF Toolbox函数有bug???
rlgc_params_even = my_s2rlgc(s_cc,LL,freq,Z0); % 不能用Zcomm!
rlgc_params_odd = my_s2rlgc(s_dd,LL,freq,Z0); % 不能用Zdiff!
gamma_e = complex(rlgc_params_even.alpha,rlgc_params_even.beta);
gamma_o = complex(rlgc_params_odd.alpha,rlgc_params_odd.beta);
Zc_e = rlgc_params_even.Zc;
Zc_o = rlgc_params_odd.Zc;

% 3. Relate the Odd and Even modes to the Differential and Common modes of
% propagation.


% 4. Using the propagation constant and the characteristic impedance for
% the Odd/Even modes, the Spice parameters are derived.
R11_cal = 0.5*real(gamma_e.*Zc_e + gamma_o.*Zc_o);
R12_cal = 0.5*real(gamma_e.*Zc_e - gamma_o.*Zc_o);
L11_cal = 0.5*imag(gamma_e.*Zc_e + gamma_o.*Zc_o)./(2*pi*freq);
L12_cal = 0.5*imag(gamma_e.*Zc_e - gamma_o.*Zc_o)./(2*pi*freq);
G11_cal = 0.5*real(gamma_e./Zc_e + gamma_o./Zc_o);
G12_cal = 0.5*real(gamma_e./Zc_e - gamma_o./Zc_o); % +-? +!
C11_cal = 0.5*imag(gamma_e./Zc_e + gamma_o./Zc_o)./(2*pi*freq);
C12_cal = 0.5*imag(gamma_e./Zc_e - gamma_o./Zc_o)./(2*pi*freq); % +-? +!

% build the  N-by-N-by-M RLCG Matrix
R_cal(2,2,:) = R11_cal;R_cal(1,1,:) = R11_cal;
R_cal(1,2,:) = R12_cal;R_cal(2,1,:) = R12_cal;
L_cal(2,2,:) = L11_cal;L_cal(1,1,:) = L11_cal;
L_cal(1,2,:) = L12_cal;L_cal(2,1,:) = L12_cal;
C_cal(2,2,:) = C11_cal;C_cal(1,1,:) = C11_cal;
C_cal(1,2,:) = C12_cal;C_cal(2,1,:) = C12_cal;
G_cal(2,2,:) = G11_cal;G_cal(1,1,:) = G11_cal;
G_cal(1,2,:) = G12_cal;G_cal(2,1,:) = G12_cal;

% Self parameters
Rs = R11_cal;
Rm = R12_cal;
Ls = L11_cal;
Lm = L12_cal;
Gs = G11_cal + G12_cal;
Gm = -G12_cal;
Cs = C11_cal + C12_cal;
Cm = -C12_cal;

%% Verify the extracted model 
s_params_cal = snp2smp(rlgc2s(R_cal,L_cal,G_cal,C_cal,LL,freq,Z0),Z0,[1 3 2 4]);
%% Verify extracted RLGC
% RLGC: even mode
figure('Name','RLGC: even mode')
subplot(221)
plot(freq,rlgc_params_even.R,'-')
hold on
plot(freq,data_SI_EM{:,6},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R(Ohm/m)');
title('R:even mode');
legend({'extracted','simulated'},'Location','best')

subplot(222)
plot(freq,rlgc_params_even.L,'-')
hold on
plot(freq,data_SI_EM{:,5},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L(H/m)');
title('L:even mode');
legend({'extracted','simulated'},'Location','best')

subplot(223)
plot(freq,rlgc_params_even.G,'-')
hold on
plot(freq,data_SI_EM{:,8},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G(S/m)');
title('G:even mode');
legend({'extracted','simulated'},'Location','best')

subplot(224)
plot(freq,rlgc_params_even.C,'-')
hold on
plot(freq,data_SI_OM{:,7},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C(F/m)');
title('C:even mode');
legend({'extracted','simulated'},'Location','best')

% RLGC: odd mode
figure('Name','RLGC: odd mode')
subplot(221)
plot(freq,rlgc_params_odd.R,'-')
hold on
plot(freq,data_SI_OM{:,6},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R(Ohm/m)');
title('R:odd mode');
legend({'extracted','simulated'},'Location','best')

subplot(222)
plot(freq,rlgc_params_odd.L,'-')
hold on
plot(freq,data_SI_OM{:,5},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L(H/m)');
title('L:odd mode');
legend({'extracted','simulated'},'Location','best')

subplot(223)
plot(freq,rlgc_params_odd.G,'-')
hold on
plot(freq,data_SI_OM{:,8},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G(S/m)');
title('G:odd mode');
legend({'extracted','simulated'},'Location','best')

subplot(224)
plot(freq,rlgc_params_odd.C,'-')
hold on
plot(freq,data_SI_OM{:,7},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C(F/m)');
title('C:odd mode');
legend({'extracted','simulated'},'Location','best')

%% Verify extracted propagation constant
% propagation constant
figure('Name','propagation constant: even/odd mode')
subplot(221)
plot(freq,real(gamma_e),'-')
hold on
plot(freq,data_SI_EM{:,16},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('\alpha(Np/m)');
title('\alpha:even mode');
legend({'extracted','simulated'},'Location','best')

subplot(222)
plot(freq,imag(gamma_e),'g-')
hold on
plot(freq,data_SI_EM{:,17},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('\beta(rad/m)');
title('\beta:even mode');
legend({'extracted','simulated'},'Location','best')

subplot(223)
plot(freq,real(gamma_o),'-')
hold on
plot(freq,data_SI_OM{:,16},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('\alpha(Np/m)');
title('\alpha:odd mode');
legend({'extracted','simulated'},'Location','best')

subplot(224)
plot(freq,imag(gamma_o),'-')
hold on
plot(freq,data_SI_OM{:,17},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('\beta(rad/m)');
title('\beta:odd mode');
legend({'extracted','simulated'},'Location','best')

%% Verify characteristic line impedances
% characteristic line impedances
figure('Name','characteristic line impedances: even/odd mode')
subplot(221)
plot(freq,real(Zc_e),'r-')
hold on
plot(freq,data_SI_EM{:,2},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('\Re(Zc,e)(Ohms)');
title('Zc:even mode');
legend({'extracted','simulated'},'Location','best')

subplot(222)
plot(freq,imag(Zc_e),'r-')
hold on
plot(freq,data_SI_EM{:,3},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('\Im(Zc,e)(Ohms)');
title('Zc:even mode');
legend({'extracted','simulated'},'Location','best')

subplot(223)
plot(freq,real(Zc_o),'r-')
hold on
plot(freq,data_SI_OM{:,2},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('\Re(Zc,o)(Ohms)');
title('Zc:odd mode');
legend({'extracted','simulated'},'Location','best')

subplot(224)
plot(freq,imag(Zc_o),'r-')
hold on
plot(freq,data_SI_OM{:,3},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('\Im(Zc,o)(Ohms)');
title('Zc:odd mode');
legend({'extracted','simulated'},'Location','best')

%% Verify spice-Parameters:tabular RLGC model(W-element)
figure('Name','spice-Parameters: R & L matrix')
subplot(221)
plot(freq,R11_cal,'g-')
hold on
plot(freq,R11_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R11(Ohms/m)');
title('R matrix');
legend({'extracted','simulated'},'Location','best')

subplot(222)
plot(freq,R12_cal,'g-')
hold on
plot(freq,R12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R12(Ohms/m)');
title('R matrix');
legend({'extracted','simulated'},'Location','best')

subplot(223)
plot(freq,L11_cal,'g-')
hold on
plot(freq,L11_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L11(H/m)');
title('L matrix');
legend({'extracted','simulated'},'Location','best')

subplot(224)
plot(freq,L12_cal,'g-')
hold on
plot(freq,L12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L12(H/m)');
title('L matrix');
legend({'extracted','simulated'},'Location','best')

figure('Name','spice-Parameters: G & C matrix')
subplot(221)
plot(freq,G11_cal,'g-')
hold on
plot(freq,G11_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G11(S/m)');
title('G matrix');
legend({'extracted','simulated'},'Location','best')

subplot(222)
plot(freq,G12_cal,'g-')
hold on
plot(freq,G12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G12(S/m)');
title('G matrix');
legend({'extracted','simulated'},'Location','best')

subplot(223)
plot(freq,C11_cal,'g-')
hold on
plot(freq,C11_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C11(F/m)');
% ylim([0 inf])
title('C matrix');
legend({'extracted','simulated'},'Location','best')

subplot(224)
plot(freq,C12_cal,'g-')
hold on
plot(freq,C12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C12(F/m)');
% ylim([-inf 0])
title('C matrix');
legend({'extracted','simulated'},'Location','best')

%% Verify Self parameters?
figure('Name','Self parameters: R & L')
subplot(221)
plot(freq,Rs,'g-')
hold on
plot(freq,R11_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('Rs(Ohms/m)');
title('Rs');
legend({'extracted','simulated'},'Location','best')

subplot(222)
plot(freq,Rm,'g-')
hold on
plot(freq,R12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('Rm(Ohms/m)');
title('Rm');
legend({'extracted','simulated'},'Location','best')

subplot(223)
plot(freq,Ls,'g-')
hold on
plot(freq,L11_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('Ls(H/m)');
title('Ls');
legend({'extracted','simulated'},'Location','best')

subplot(224)
plot(freq,Lm,'g-')
hold on
plot(freq,L12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('Lm(H/m)');
title('Lm');
legend({'extracted','simulated'},'Location','best')

figure('Name','Self parameters: G & C')
subplot(221)
plot(freq,Gs,'g-')
hold on
plot(freq,G11_sim + G12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('Gs(S/m)');
title('Gs');
legend({'extracted','simulated'},'Location','best')

subplot(222)
plot(freq,Gm,'g-')
hold on
plot(freq,-G12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('Gm(S/m)');
title('Gm');
legend({'extracted','simulated'},'Location','best')

subplot(223)
plot(freq,Cs,'g-')
hold on
plot(freq,C11_sim+C12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('Cs(F/m)');
title('Cs');
legend({'extracted','simulated'},'Location','best')

subplot(224)
plot(freq,Cm,'g-')
hold on
plot(freq,-C12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('Cm(F/m)');
title('Cm');
legend({'extracted','simulated'},'Location','best')

%% Test: MATLAB function "s2rlgc" performed at 4-port S directly
figure('Name','"s2rlgc(Matlab)", spice-Parameters: R & L matrix')
subplot(221)
plot(freq,squeeze(rlgc_params_matlab.R(1,1,:)),'g-')
hold on
plot(freq,R11_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R11(Ohms/m)');
title('R matrix');
legend({'s2rlgc(Matlab)','simulated'},'Location','best')

subplot(222)
plot(freq,squeeze(rlgc_params_matlab.R(1,2,:)),'g-')
hold on
plot(freq,R12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('R12(Ohms/m)');
title('R matrix');
legend({'s2rlgc(Matlab)','simulated'},'Location','best')

subplot(223)
plot(freq,squeeze(rlgc_params_matlab.L(1,1,:)),'g-')
hold on
plot(freq,L11_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L11(H/m)');
title('L matrix');
legend({'s2rlgc(Matlab)','simulated'},'Location','best')

subplot(224)
plot(freq,squeeze(rlgc_params_matlab.L(1,2,:)),'g-')
hold on
plot(freq,L12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('L12(H/m)');
title('L matrix');
legend({'s2rlgc(Matlab)','simulated'},'Location','best')

figure('Name','"s2rlgc(Matlab)", spice-Parameters: G & C matrix')
subplot(221)
plot(freq,squeeze(rlgc_params_matlab.G(1,1,:)),'g-')
hold on
plot(freq,G11_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G11(S/m)');
title('G matrix');
legend({'s2rlgc(Matlab)','simulated'},'Location','best')

subplot(222)
plot(freq,squeeze(rlgc_params_matlab.G(1,2,:)),'g-')
hold on
plot(freq,G12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('G12(S/m)');
title('G matrix');
legend({'s2rlgc(Matlab)','simulated'},'Location','best')

subplot(223)
plot(freq,squeeze(rlgc_params_matlab.C(1,1,:)),'g-')
hold on
plot(freq,C11_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C11(F/m)');
title('C matrix');
legend({'s2rlgc(Matlab)','simulated'},'Location','best')

subplot(224)
plot(freq,squeeze(rlgc_params_matlab.C(1,2,:)),'g-')
hold on
plot(freq,C12_sim,'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('C12(F/m)');
title('C matrix');
legend({'s2rlgc(Matlab)','simulated'},'Location','best')

%% Test: MATLAB function "rlgc2s" performed after "s2rlgc"
figure('Name','"s2rlgc2s(Matlab)", S parameters')
subplot(221)
plot(freq,db(squeeze(s_params_matlab(1,1,:))),'g-')
hold on
plot(freq,db(squeeze(s_params(1,1,:))),'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('S11(dB)');
title('S11');
legend({'s2rlgc2s(Matlab)','simulated'},'Location','best')

subplot(222)
plot(freq,db(squeeze(s_params_matlab(1,2,:))),'g-')
hold on
plot(freq,db(squeeze(s_params(1,2,:))),'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('S12(dB)');
title('S11');
legend({'s2rlgc2s(Matlab)','simulated'},'Location','best')

subplot(223)
plot(freq,db(squeeze(s_params_matlab(1,3,:))),'g-')
hold on
plot(freq,db(squeeze(s_params(1,3,:))),'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('S13(dB)');
title('S13');
legend({'s2rlgc2s(Matlab)','simulated'},'Location','best')

subplot(224)
plot(freq,db(squeeze(s_params_matlab(1,4,:))),'g-')
hold on
plot(freq,db(squeeze(s_params(1,4,:))),'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('S14(dB)');
title('S14');
legend({'s2rlgc2s(Matlab)','simulated'},'Location','best')

%% Verify extracted model: Using MATLAB function "rlgc2s"
figure('Name','"Calculated [S] using extracted RLGC model"')
subplot(221)
plot(freq,db(squeeze(s_params_cal(1,1,:))),'g-')
hold on
plot(freq,db(squeeze(s_params(1,1,:))),'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('S11(dB)');
title('S11');
legend({'extracted','simulated'},'Location','best')

subplot(222)
plot(freq,db(squeeze(s_params_cal(1,2,:))),'g-')
hold on
plot(freq,db(squeeze(s_params(1,2,:))),'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('S12(dB)');
title('S12');

subplot(223)
plot(freq,db(squeeze(s_params_cal(1,3,:))),'g-')
hold on
plot(freq,db(squeeze(s_params(1,3,:))),'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('S13(dB)');
title('S13');
legend({'extracted','simulated'},'Location','best')

subplot(224)
plot(freq,db(squeeze(s_params_cal(1,4,:))),'g-')
hold on
plot(freq,db(squeeze(s_params(1,4,:))),'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('S14(dB)');
title('S14');
legend({'extracted','simulated'},'Location','best')

%% Test: prop const calc by RF Toolbox function "s2rlgc" 
figure('Name','propagation constant: 4 Port')
subplot(221)
plot(freq,real(gamma_e),'-')
hold on
plot(freq,data_SI_EM{:,16},'k--')
hold off
grid on
xlabel('Freq(Hz)');
ylabel('\alpha(Np/m)');
title('\alpha:even mode');
legend({'extracted','simulated'},'Location','best')

%% Test: Consistent of function "s2rlgc"
% 
