
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