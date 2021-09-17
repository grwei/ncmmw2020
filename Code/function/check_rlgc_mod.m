function [num_lines, freq_pts] = check_rlgc_mod(allR,allL,allG,allC,freq)
% Test the properties of RLGC matrices 

% Must be numeric
if (~isnumeric(allR) || isempty(allR) || ...
        any(any(any(isnan(allR)))) || any(any(any(isinf(allR)))))
    error(message('rf:check_rlgc:isnumericR'));
end

if (~isnumeric(allL) || isempty(allL) || ...
        any(any(any(isnan(allL)))) || any(any(any(isinf(allL)))))
    error(message('rf:check_rlgc:isnumericL'));
end

if (~isnumeric(allG) || isempty(allG) || ...
        any(any(any(isnan(allG)))) || any(any(any(isinf(allG)))))
    error(message('rf:check_rlgc:isnumericG'));
end

if (~isnumeric(allC) || isempty(allC) || ...
        any(any(any(isnan(allC)))) || any(any(any(isinf(allC)))))
    error(message('rf:check_rlgc:isnumericC'));
end

if((size(allR) == size(allL)) & (size(allG) == size(allC)) ...
        & (size(allR) == size(allC))) %#ok<AND2>
    num_lines = size(allR,1);       % Number of transmission lines
    freq_pts  = size(allR,3);       % number of frequency points
else
    error(message('rf:check_rlgc:samesizeRLGC'));
end

num_terms = num_lines*num_lines;
for m=1:freq_pts    
    R = allR(:,:,m);
    L = allL(:,:,m);
    G = allG(:,:,m);
    C = allC(:,:,m);
    
    % All matrices are symmetric. 
    if( ~isequal(R, R.') || ~isequal(L, L.') || ~isequal(G, G.') || ~isequal(C, C.'))
        error(message('rf:check_rlgc:symmetric'));
    end
    
    % The diagonal terms of L and C are positive, non-zero. 
    if((any(diag(L)<= 0) ||  any(diag(C)<= 0)) &&  freq(m) ~= 0.0)
        warning(message('rf:check_rlgc:diagLC'));  % mod: error -> warnning, Wei, 20191023
    end
    
    % The diagonal terms of R and G are non-negative (can be zero). 
    if(any(diag(R)< 0) ||  any(diag(G)< 0))
        warning(message('rf:check_rlgc:diagRG'));  % mod: error -> warnning, Wei, 20191019        
    end
        
    % Off-diagonal terms of the L matrix are non-negative.
    if(any(any(L <0)) && num_lines>1)
        warning(message('rf:check_rlgc:offdiagL')); % mod: error -> warnning, Wei, 20191019
    end
        
    % Off-diagonal terms of C and G matrices are non-positive. 
    if((sum(sum(C <=0)) == num_terms*(num_terms-1)) && num_lines>1)
        warning(message('rf:check_rlgc:offdiagC')); % mod: error -> warnning, Wei, 20191019         
    end
    
    if((sum(sum(G <=0)) == num_terms*(num_terms-1)) && num_lines>1)
        error(message('rf:check_rlgc:offdiagG'));
    end
end

end