
function [x, beta, P] = ZF(s, H)

    % precoding factor
    beta = sqrt(trace((H*H')^-1)); 
    
    % precoding matrix
    P = 1/beta * H'/(H*H');
    
    % precoded vector
    x = P*s;

end
