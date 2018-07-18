We borrow this code from https://github.com/quantizedmassivemimo/1bit_precoding;

function [x, beta, P] = MRT(s, H)

    % number BS antennas
    [~, B] = size(H);

    % precoding factor
    beta = sqrt(trace(H'*H))/B;
    
    % precoding matrix
    P = 1/B/beta * H';
    
    % precoded vector
    x = P*s;
                
end
