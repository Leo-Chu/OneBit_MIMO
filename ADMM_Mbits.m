
function [x, beta,xRest] = ADMM_Mbits(par,s,H,N0)

    HR0 = [ real(H) -imag(H) ; imag(H) real(H) ];
    sR = [real(s) ; imag(s)];
    
    if par.b == 1
        C = eye(2*par.N);  bps = 2; 
    elseif par.b == 2
        Q1 = sqrt(5);    
%         Q = 20*log10(Q1);
        C = [2*eye(2*par.N) eye(2*par.N)]/Q1; bps = 4; %N0 = N0*Q1; % ((bps/2)^2+1)
    elseif par.b == 3
        Q1 = sqrt(21);
%         Q = 20*log10(Q1);
        C = [4*eye(2*par.N) 2*eye(2*par.N) eye(2*par.N)]/Q1; bps = 6; %N0 = N0*Q1;
    else
        disp('Not supportted at current version!!!');
    end
    HR = HR0*C;

    x = zeros(par.N*bps,1);
    y = zeros(par.N*bps,1);
    
    if(par.N<64)
          gain = .05;
    else
          gain = 1;  
    end
%     gain = 1/sqrt(par.U*par.N);
    epsilon = 1e-5; 
    
    A = HR'*HR + 0.5/gain*eye(par.N*bps);
    sREG = A\(HR'*sR);
    
    iter = 1e2;
    for t=1:iter
        u = sREG + 0.5/gain*(A\(2*x-y));
        xold = x;
        x = prox_inf_norm(y+u-x,4*par.U*par.N*N0);        
        if norm(x-xold)/norm(x)<epsilon
            break;
        end
        y = y +  u- x;
    end
%     disp(t)
    
    xRest = sign(x);
%     if par.b == 1
%           Q1 = 1;    %Q = Q1;  
%     elseif par.b == 2
%           Q1 = sqrt(5);  %Q = 20*log10(Q1);
%     else
%           Q1 = sqrt(21); %Q = 20*log10(Q1);
%     end
%     x0 = C*xRest/Q1; 
    
    x0 = C*xRest;     
    x = 1/sqrt(par.N*2)*(x0(1:par.N,1)+1i*x0(par.N+1:2*par.N,1)); %*Q  *Q*par.N 2*
    
    beta = real(x'*H'*s)/(norm(H*x,2)^2+par.U*N0);   
   
    if beta < 0
        x = -x;
        beta = -beta;
    end

end

%% revised from http://standford.edu/~boyd/papers/prox_algs.html

function xk = prox_inf_norm(u,lambda)

        N = length(u);
        u_abs = abs(u);
        us = (cumsum(sort(u_abs,'descend')))./(lambda+(1:N)');
        alphaopt = max(us);
%         xk = min(u_abs,alphaopt).*sign(u); 
        if alphaopt>0 
                xk = min(u_abs,alphaopt).*sign(u); % truncation step
        else
                xk = zeros(size(u)); % if t is big, then solution is zero
        end  
end

