function [output] = zGamma2rlgc(Zc,gamma,freq)
%ZGAMMA2RLGC convert complex propagation constant and characteristic
%impedance to RLGC
% 
%   Zc is a complex N-by-N-by-M Characteristic line impedance(ohm)
%   gamma = complex(alpha,beta),where
%   alpha is a real N-by-N-by-M attenuation constant (Nepers/m)
%   beta is a real N-by-N-by-M phase constant (radians/m)
%   FREQ is a real Mx1 frequency vector
% 
%   The outputs are per unit length transmission line parameters
%   OUTPUT.R is a real N-by-N-by-M Resistance matrix (ohm/m)
%   OUTPUT.L is a real N-by-N-by-M Inductance matrix (H/m)
%   OUTPUT.C is a real N-by-N-by-M Capacitance matrix (F/m)
%   OUTPUT.G is a real N-by-N-by-M Conductance matrix (S/m)

%% Allocate the memory for RLGC matrix
numLines = size(gamma,1);
freqpts = size(Zc,3);
R   = zeros(numLines,numLines,freqpts);      % Resistance matrix (ohm/m)
L   = R;                                     % Inductance matrix (H/m)
G   = R;                                     % Conductance matrix (S/m)
C   = R;                                     % Capacitance matrix (F/m)

%% Calculate the RLGC matrices
for m= 1:freqpts
    R_temp = real(Zc(:,:,m)*gamma(:,:,m));
    L_temp = imag(Zc(:,:,m)*gamma(:,:,m))./(2*pi*freq(m));
    G_temp = real(gamma(:,:,m)/Zc(:,:,m));
    C_temp = imag(gamma(:,:,m)/Zc(:,:,m))./(2*pi*freq(m));
    
    R(:,:,m) = 0.5*(R_temp + R_temp.');
    L(:,:,m) = 0.5*(L_temp + L_temp.');
    G(:,:,m) = 0.5*(G_temp + G_temp.');
    C(:,:,m) = 0.5*(C_temp + C_temp.');
end
output = struct('R',R,'L',L,'G',G,'C',C);

end
