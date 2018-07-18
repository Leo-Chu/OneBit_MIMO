We borrow this code from https://github.com/quantizedmassivemimo/1bit_precoding;

function [x, beta, P] = WF(s, H, N0)

    % number of UEs
    [U, ~] = size(H);
    
    % precoding matrix (before normalization)
    T = H' / (H*H' + U*N0*eye(U));

    % precoding factor
    beta = sqrt(real(trace((T*T'))));
    
    % precoding matrix
    P = 1/beta*T;
    
    % precoded vector
    x = P*s;

end
