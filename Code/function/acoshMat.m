function [acoshMat] = acoshMat(mat)
%ACOSHMAT inverse hyperbolic cosine of a N¡ÁN¡ÁM matrix
%   reference: Matlab RF Toolbox function: s2rlgc.m

numLines = size(mat,1);
freqpts = size(mat,3);

% Calculate the eigenvalues and eigenvectors for Td
V = zeros(numLines,numLines,freqpts);      % Right Eigenvectors of Td
D = zeros(numLines,freqpts);                 % Eigenvalues of Td
for idx=1:freqpts
    [V(:,:,idx),D(:,idx)] = eig(mat(:,:,idx),'vector');
end

% Method 2
% Compute quantities with diagonalization from sort
gammaTrLen = acosh(D);
% Eigen values need to be sorted at each frequency
[~, idxGammaTrLenSorted] = sort(abs(real(gammaTrLen)));
gammaTrLenSorted = cell2mat(arrayfun(                               ...
    @(idx) gammaTrLen(idxGammaTrLenSorted(:,idx),idx), 1:freqpts,   ...
    'UniformOutput', false));
gammaTrLenSortedUnwrap = real(gammaTrLenSorted)    +                ...
    1i*unwrap(imag(gammaTrLenSorted),[],2);
Vsorted = reshape(cell2mat(                                         ...
    arrayfun(@(idx) V(:,idxGammaTrLenSorted(:,idx),idx),1:freqpts,  ...
    'UniformOutput', false)),numLines,numLines,freqpts);
gamma2 = zeros(numLines,numLines,freqpts);
for idx=1:freqpts
    gamma2(:,:,idx) = Vsorted(:,:,idx) * ...
        (diag(gammaTrLenSortedUnwrap(:,idx)))/Vsorted(:,:,idx);
end

alpha = real(gamma2);
beta  = imag(gamma2);

for m = 1:numLines
    for n = 1:numLines
        beta(m,n,:)  = abs(beta(m,n,:));
        alpha(m,n,:) = abs(alpha(m,n,:));
    end
end
acoshMat = complex(alpha,beta);

end

