function [rlgc_params] = mms2rlgc(smm,LL,freq,Z0)
%Extracting the RLCG parameters from MMS for differential TL
%   [rlgc_params] = mms2rlgc(smm,LL,freq,Z0)
%   The outputs are per unit length transmission line parameters
%   OUTPUT.R is a real N-by-N-by-M Resistance matrix (ohm/m)
%   OUTPUT.L is a real N-by-N-by-M Inductance matrix (H/m)
%   OUTPUT.C is a real N-by-N-by-M Capacitance matrix (F/m)
%   OUTPUT.G is a real N-by-N-by-M Conductance matrix (S/m)

%% Calculate the propagation constant and the characteristic impedance
rlgc_params_even = s2rlgc_mod2(smm.cc,LL,freq,Z0); % 不能用Zcomm!
rlgc_params_odd = s2rlgc_mod2(smm.dd,LL,freq,Z0); % 不能用Zdiff!
gamma_e = complex(rlgc_params_even.alpha,rlgc_params_even.beta);
gamma_o = complex(rlgc_params_odd.alpha,rlgc_params_odd.beta);
Zc_e = rlgc_params_even.Zc;
Zc_o = rlgc_params_odd.Zc;

%% Using the propagation constant and the characteristic impedance for
% the Odd/Even modes, the Spice parameters are derived.
R11 = 0.5*real(gamma_e.*Zc_e + gamma_o.*Zc_o);
R12 = 0.5*real(gamma_e.*Zc_e - gamma_o.*Zc_o);
L11 = 0.5*imag(gamma_e.*Zc_e + gamma_o.*Zc_o)./(2*pi*freq);
L12 = 0.5*imag(gamma_e.*Zc_e - gamma_o.*Zc_o)./(2*pi*freq);
G11 = 0.5*real(gamma_e./Zc_e + gamma_o./Zc_o);
G12 = 0.5*real(gamma_e./Zc_e - gamma_o./Zc_o); 
C11 = 0.5*imag(gamma_e./Zc_e + gamma_o./Zc_o)./(2*pi*freq);
C12 = 0.5*imag(gamma_e./Zc_e - gamma_o./Zc_o)./(2*pi*freq); 

% build the  N-by-N-by-M RLCG Matrix
R(2,2,:) = R11;R(1,1,:) = R11;
R(1,2,:) = R12;R(2,1,:) = R12;
L(2,2,:) = L11;L(1,1,:) = L11;
L(1,2,:) = L12;L(2,1,:) = L12;
C(2,2,:) = C11;C(1,1,:) = C11;
C(1,2,:) = C12;C(2,1,:) = C12;
G(2,2,:) = G11;G(1,1,:) = G11;
G(1,2,:) = G12;G(2,1,:) = G12;

% struct
rlgc_params.R = R;
rlgc_params.L = L;
rlgc_params.G = G;
rlgc_params.C = C;
rlgc_params.even = rlgc_params_even;
rlgc_params.odd = rlgc_params_odd;

%% debug(Remark this section after debugging)
debug = false;
if debug == true
    %% Verify extracted propagation constant
    % propagation constant
    figure('Name','Debug: propagation constant')
    subplot(221)
    plot(freq,real(gamma_e))
    grid on
    xlabel('Freq(Hz)');
    ylabel('\alpha(Np/m)');
    title('\alpha:even mode');
    legend({'\alpha_e'},'Location','best')
    
    subplot(222)
    plot(freq,imag(gamma_e))
    grid on
    xlabel('Freq(Hz)');
    ylabel('\beta(Rad/m)');
    title('\beta:even mode');
    legend({'\beta_e'},'Location','best')
    
    subplot(223)
    plot(freq,real(gamma_o))
    grid on
    xlabel('Freq(Hz)');
    ylabel('\alpha(Np/m)');
    title('\alpha:odd mode');
    legend({'\alpha_o'},'Location','best')
    
    subplot(224)
    plot(freq,imag(gamma_o))
    grid on
    xlabel('Freq(Hz)');
    ylabel('\beta(Rad/m)');
    title('\beta:odd mode');
    legend({'\beta_o'},'Location','best')
end

end

