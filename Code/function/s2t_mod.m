function [t_params] = s2t_mod(s_params,z0)
%S2T_MOD convert S matrix to T matrix
%by s -> z -> t, different from the direct conversion in function"s2abcd"
%   [t_params] = s2t_mod(s_params,z0)

numLines = size(s_params,1)/2;
%% Convert S matrix to Z matrix
% Allocate memory for the Z-parameters
z_params = zeros(size(s_params)); 
freqpts = size(z_params,3);
% Calc the Z-parameters: Z = Z0 * (I + S) * inv(I - S)
I = eye(size(s_params, 1));
for k = 1:freqpts
    z_params(:,:,k) = (z0 * (I + s_params(:,:,k))) /(I - s_params(:,:,k));
end

%% Convert Z matrix to T matrix
Z11 = z_params(1:numLines,1:numLines,:);
Z12 = z_params(1:numLines,numLines+1:end,:);
Z21 = z_params(numLines+1:end,1:numLines,:);
Z22 = z_params(numLines+1:end,numLines+1:end,:);

TD = zeros(size(Z11));
TA = TD;
TB = TD;
TC = TD;
t_params = zeros(size(s_params));
for k = 1:freqpts
    TA(:,:,k) = Z11(:,:,k)/Z21(:,:,k);
    TB(:,:,k) = Z11(:,:,k)/Z21(:,:,k)*Z22(:,:,k)-Z12(:,:,k);
    TC(:,:,k) = Z21(:,:,k)\eye(numLines);
    TD(:,:,k) = Z21(:,:,k)\Z22(:,:,k);
    t_params(:,:,k) = [TA(:,:,k),TB(:,:,k);TC(:,:,k),TD(:,:,k)];
end

end

