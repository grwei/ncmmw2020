function [rlgc_params_fit] = rlgc_fit(rlgc_params,freq,outliers)
%RLGC_FIT Fit the rlgc params
%   [rlgc_params_fit] = rlgc_fit(rlgc_params,freq,outliers)
%% Fitting the RLCG parameters obtained by the closed-form equations with the model
% R matrix
R = zeros(size(rlgc_params.R));
[R11_fitresult, R11_gof] = R_fit(freq', squeeze(rlgc_params.R(1,1,:))', outliers);
R(1,1,:) = R11_fitresult(freq');
[R12_fitresult, R12_gof] = R_fit(freq', squeeze(rlgc_params.R(1,2,:))', outliers);
R(1,2,:) = R12_fitresult(freq');
[R13_fitresult, R13_gof] = R_fit(freq', squeeze(rlgc_params.R(1,3,:))', outliers);
R(1,3,:) = R13_fitresult(freq');
[R14_fitresult, R14_gof] = R_fit(freq', squeeze(rlgc_params.R(1,4,:))', outliers);
R(1,4,:) = R14_fitresult(freq');
[R21_fitresult, R21_gof] = R_fit(freq', squeeze(rlgc_params.R(2,1,:))', outliers);
R(2,1,:) = R21_fitresult(freq');
[R22_fitresult, R22_gof] = R_fit(freq', squeeze(rlgc_params.R(2,2,:))', outliers);
R(2,2,:) = R22_fitresult(freq');
[R23_fitresult, R23_gof] = R_fit(freq', squeeze(rlgc_params.R(2,3,:))', outliers);
R(2,3,:) = R23_fitresult(freq');
[R24_fitresult, R24_gof] = R_fit(freq', squeeze(rlgc_params.R(2,4,:))', outliers);
R(2,4,:) = R24_fitresult(freq');
[R31_fitresult, R31_gof] = R_fit(freq', squeeze(rlgc_params.R(3,1,:))', outliers);
R(3,1,:) = R31_fitresult(freq');
[R32_fitresult, R32_gof] = R_fit(freq', squeeze(rlgc_params.R(3,2,:))', outliers);
R(3,2,:) = R32_fitresult(freq');
[R33_fitresult, R33_gof] = R_fit(freq', squeeze(rlgc_params.R(3,3,:))', outliers);
R(3,3,:) = R33_fitresult(freq');
[R34_fitresult, R34_gof] = R_fit(freq', squeeze(rlgc_params.R(3,4,:))', outliers);
R(3,4,:) = R34_fitresult(freq');
[R41_fitresult, R41_gof] = R_fit(freq', squeeze(rlgc_params.R(4,1,:))', outliers);
R(4,1,:)= R41_fitresult(freq');
[R42_fitresult, R42_gof] = R_fit(freq', squeeze(rlgc_params.R(4,2,:))', outliers);
R(4,2,:) = R42_fitresult(freq');
[R43_fitresult, R43_gof] = R_fit(freq', squeeze(rlgc_params.R(4,3,:))', outliers);
R(4,3,:) = R43_fitresult(freq');
[R44_fitresult, R44_gof] = R_fit(freq', squeeze(rlgc_params.R(4,4,:))', outliers);
R(4,4,:) = R44_fitresult(freq');

% L matrix
L = zeros(size(rlgc_params.L));
[L11_fitresult, L11_gof] = L_fit(freq', squeeze(rlgc_params.L(1,1,:))', outliers);
L(1,1,:) = L11_fitresult(freq');
[L12_fitresult, L12_gof] = L_fit(freq', squeeze(rlgc_params.L(1,2,:))', outliers);
L(1,2,:) = L12_fitresult(freq');
[L13_fitresult, L13_gof] = L_fit(freq', squeeze(rlgc_params.L(1,3,:))', outliers);
L(1,3,:) = L13_fitresult(freq');
[L14_fitresult, L14_gof] = L_fit(freq', squeeze(rlgc_params.L(1,4,:))', outliers);
L(1,4,:) = L14_fitresult(freq');
[L21_fitresult, L21_gof] = L_fit(freq', squeeze(rlgc_params.L(2,1,:))', outliers);
L(2,1,:) = L21_fitresult(freq');
[L22_fitresult, L22_gof] = L_fit(freq', squeeze(rlgc_params.L(2,2,:))', outliers);
L(2,2,:) = L22_fitresult(freq');
[L23_fitresult, L23_gof] = L_fit(freq', squeeze(rlgc_params.L(2,3,:))', outliers);
L(2,3,:) = L23_fitresult(freq');
[L24_fitresult, L24_gof] = L_fit(freq', squeeze(rlgc_params.L(2,4,:))', outliers);
L(2,4,:) = L24_fitresult(freq');
[L31_fitresult, L31_gof] = L_fit(freq', squeeze(rlgc_params.L(3,1,:))', outliers);
L(3,1,:) = L31_fitresult(freq');
[L32_fitresult, L32_gof] = L_fit(freq', squeeze(rlgc_params.L(3,2,:))', outliers);
L(3,2,:) = L32_fitresult(freq');
[L33_fitresult, L33_gof] = L_fit(freq', squeeze(rlgc_params.L(3,3,:))', outliers);
L(3,3,:) = L33_fitresult(freq');
[L34_fitresult, L34_gof] = L_fit(freq', squeeze(rlgc_params.L(3,4,:))', outliers);
L(3,4,:) = L34_fitresult(freq');
[L41_fitresult, L41_gof] = L_fit(freq', squeeze(rlgc_params.L(4,1,:))', outliers);
L(4,1,:)= L41_fitresult(freq');
[L42_fitresult, L42_gof] = L_fit(freq', squeeze(rlgc_params.L(4,2,:))', outliers);
L(4,2,:) = L42_fitresult(freq');
[L43_fitresult, L43_gof] = L_fit(freq', squeeze(rlgc_params.L(4,3,:))', outliers);
L(4,3,:) = L43_fitresult(freq');
[L44_fitresult, L44_gof] = L_fit(freq', squeeze(rlgc_params.L(4,4,:))', outliers);
L(4,4,:) = L44_fitresult(freq');

% C matrix
C = zeros(size(rlgc_params.C));
[C11_fitresult, C11_gof] = C_fit(freq', squeeze(rlgc_params.C(1,1,:))', outliers);
C(1,1,:) = C11_fitresult(freq');
[C12_fitresult, C12_gof] = C_fit(freq', squeeze(rlgc_params.C(1,2,:))', outliers);
C(1,2,:) = C12_fitresult(freq');
[C13_fitresult, C13_gof] = C_fit(freq', squeeze(rlgc_params.C(1,3,:))', outliers);
C(1,3,:) = C13_fitresult(freq');
[C14_fitresult, C14_gof] = C_fit(freq', squeeze(rlgc_params.C(1,4,:))', outliers);
C(1,4,:) = C14_fitresult(freq');
[C21_fitresult, C21_gof] = C_fit(freq', squeeze(rlgc_params.C(2,1,:))', outliers);
C(2,1,:) = C21_fitresult(freq');
[C22_fitresult, C22_gof] = C_fit(freq', squeeze(rlgc_params.C(2,2,:))', outliers);
C(2,2,:) = C22_fitresult(freq');
[C23_fitresult, C23_gof] = C_fit(freq', squeeze(rlgc_params.C(2,3,:))', outliers);
C(2,3,:) = C23_fitresult(freq');
[C24_fitresult, C24_gof] = C_fit(freq', squeeze(rlgc_params.C(2,4,:))', outliers);
C(2,4,:) = C24_fitresult(freq');
[C31_fitresult, C31_gof] = C_fit(freq', squeeze(rlgc_params.C(3,1,:))', outliers);
C(3,1,:) = C31_fitresult(freq');
[C32_fitresult, C32_gof] = C_fit(freq', squeeze(rlgc_params.C(3,2,:))', outliers);
C(3,2,:) = C32_fitresult(freq');
[C33_fitresult, C33_gof] = C_fit(freq', squeeze(rlgc_params.C(3,3,:))', outliers);
C(3,3,:) = C33_fitresult(freq');
[C34_fitresult, C34_gof] = C_fit(freq', squeeze(rlgc_params.C(3,4,:))', outliers);
C(3,4,:) = C34_fitresult(freq');
[C41_fitresult, C41_gof] = C_fit(freq', squeeze(rlgc_params.C(4,1,:))', outliers);
C(4,1,:)= C41_fitresult(freq');
[C42_fitresult, C42_gof] = C_fit(freq', squeeze(rlgc_params.C(4,2,:))', outliers);
C(4,2,:) = C42_fitresult(freq');
[C43_fitresult, C43_gof] = C_fit(freq', squeeze(rlgc_params.C(4,3,:))', outliers);
C(4,3,:) = C43_fitresult(freq');
[C44_fitresult, C44_gof] = C_fit(freq', squeeze(rlgc_params.C(4,4,:))', outliers);
C(4,4,:) = C44_fitresult(freq');

% G matrix
G = zeros(size(rlgc_params.G));
[G11_fitresult, G11_gof] = G_fit(freq', squeeze(rlgc_params.G(1,1,:))', outliers);
G(1,1,:) = G11_fitresult(freq');
[G12_fitresult, G12_gof] = G_fit(freq', squeeze(rlgc_params.G(1,2,:))', outliers);
G(1,2,:) = G12_fitresult(freq');
[G13_fitresult, G13_gof] = G_fit(freq', squeeze(rlgc_params.G(1,3,:))', outliers);
G(1,3,:) = G13_fitresult(freq');
[G14_fitresult, G14_gof] = G_fit(freq', squeeze(rlgc_params.G(1,4,:))', outliers);
G(1,4,:) = G14_fitresult(freq');
[G21_fitresult, G21_gof] = G_fit(freq', squeeze(rlgc_params.G(2,1,:))', outliers);
G(2,1,:) = G21_fitresult(freq');
[G22_fitresult, G22_gof] = G_fit(freq', squeeze(rlgc_params.G(2,2,:))', outliers);
G(2,2,:) = G22_fitresult(freq');
[G23_fitresult, G23_gof] = G_fit(freq', squeeze(rlgc_params.G(2,3,:))', outliers);
G(2,3,:) = G23_fitresult(freq');
[G24_fitresult, G24_gof] = G_fit(freq', squeeze(rlgc_params.G(2,4,:))', outliers);
G(2,4,:) = G24_fitresult(freq');
[G31_fitresult, G31_gof] = G_fit(freq', squeeze(rlgc_params.G(3,1,:))', outliers);
G(3,1,:) = G31_fitresult(freq');
[G32_fitresult, G32_gof] = G_fit(freq', squeeze(rlgc_params.G(3,2,:))', outliers);
G(3,2,:) = G32_fitresult(freq');
[G33_fitresult, G33_gof] = G_fit(freq', squeeze(rlgc_params.G(3,3,:))', outliers);
G(3,3,:) = G33_fitresult(freq');
[G34_fitresult, G34_gof] = G_fit(freq', squeeze(rlgc_params.G(3,4,:))', outliers);
G(3,4,:) = G34_fitresult(freq');
[G41_fitresult, G41_gof] = G_fit(freq', squeeze(rlgc_params.G(4,1,:))', outliers);
G(4,1,:)= G41_fitresult(freq');
[G42_fitresult, G42_gof] = G_fit(freq', squeeze(rlgc_params.G(4,2,:))', outliers);
G(4,2,:) = G42_fitresult(freq');
[G43_fitresult, G43_gof] = G_fit(freq', squeeze(rlgc_params.G(4,3,:))', outliers);
G(4,3,:) = G43_fitresult(freq');
[G44_fitresult, G44_gof] = G_fit(freq', squeeze(rlgc_params.G(4,4,:))', outliers);
G(4,4,:) = G44_fitresult(freq');

% reciporcal (Xij=Xji).
for idx = 1:length(freq)
    R(:,:,idx) = 0.5 * (R(:,:,idx) + R(:,:,idx)'); 
    L(:,:,idx) = 0.5 * (L(:,:,idx) + L(:,:,idx)'); 
    G(:,:,idx) = 0.5 * (G(:,:,idx) + G(:,:,idx)'); 
    C(:,:,idx) = 0.5 * (C(:,:,idx) + C(:,:,idx)'); 
end

rlgc_params_fit = struct('R',R,'L',L,'G',G,'C',C);

end

