function [revGamma] = mainpha2contpha(freq,Gamma,thres)
% This function is an extension of the MATLAB built-in function unwrap,
% which Correct phase angles to produce smoother phase plots.
% [input] Gamma is a complex vector with main image part,
% thres is the tolerance of nonlinearity,range 0~1,
% revGamma is a complex vector with continuous image part,
% [output] imag(Gamma(freq==0)) should be near 0.

% Author Name: Wei Guorui
% Instructor Name: Xia Bin
% Created in: 2019-05-21 
% Last modified: 2019-05-22 04:17

imag_Gamma = imag(Gamma); % image part of Gamma
maxVal = max(imag_Gamma); 
minVal = min(imag_Gamma);
% scaled,ref to help:unwrap
lbound = -2*pi;
ubound = 2*pi;
imag_scaled = rescale(imag_Gamma,lbound,ubound); 
% unwraps radian phases by changing absolute jumps greater than pi to their 2*pi complement.
uw_scaled = unwrap(imag_scaled); 
% consider: imag(Gamma(freq==0)) should be near 0.
[p,~,mu] = polyfit(freq,uw_scaled,1); % ref to HELP POLYFIT
intercept = polyval(p,0,[],mu); % the longitudinal intercept
k = round(intercept/pi);  
% accept range: intercept/pi in range [k-tol,k+tol],where tol less than 1.
if abs(intercept/pi - k) >= thres
    warning('Nonlinearity exceeds threshold, try simulating from a lower start frequency or using a finer grid.')
end
imag_cont = (uw_scaled - k*pi - lbound)*(maxVal - minVal)/(ubound - lbound) + minVal;
revGamma = complex(real(Gamma),imag_cont);

end

