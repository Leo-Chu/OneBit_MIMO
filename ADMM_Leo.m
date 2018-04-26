%% developed by Lei Chu and Fei Wen

function [x, beta, zr,vr] = ADMM_Leo(par,s,H,N0)

    HR0 = [ real(H) -imag(H) ; imag(H) real(H) ];
    sR = [real(s) ; imag(s)];
    
    if par.b == 1
        C = eye(2*par.N);  bps = 2; %Q = bps/2; N0 = N0*Q; % 这里没有 虚拟向量扩充，所以信噪比没有变化！！
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

%         2*norm(HR'*HR,2)
    
    v = zeros(par.N*bps,1);
    z = zeros(par.N*bps,1);
    w = zeros(par.N*bps,1);
    
    iter = 5e2; c = par.U*N0; epsilon = 1e-5;

%     hv =  inv(2*HR'*HR + (2*c+rho)*eye(par.N*bps)); 
   vr=[]; zr = [];
   rho_0=1;
   rho = rho_0;
   rho_t = 2.1*norm(HR'*HR,2); 
   [uu,ss,~] = svd(HR'*HR);
   Hs = HR'*sR;
% rho = rho_t;
    for k = 1:iter
          vm1 = v;
          zm1 = z; 
          if rho<rho_t
                rho = rho_0*1.15^k;
          end
          gg = 2*diag(ss) + (2*c+rho);
%           v = uu*(diag(1./gg)*(uu'*(2*Hs + rho*z + w)));
          v = uu*((uu'*(2*Hs + rho*z + w))./gg);

          %v =  inv(2*HR'*HR + (2*c+rho)*eye(par.N*bps))*(2*HR'*sR + rho*z + w);  
          %v =  (2*HR'*HR + (2*c+rho)*eye(par.N*bps))\(2*HR'*sR + rho*z + w);  

          z = sign(v - w/rho)*norm(v - w/rho,1)/(par.N*bps);
          w = w - rho*(v - z);
          
          vr = [vr, norm(v-vm1)];
          zr = [zr, norm(z-zm1)];
%           if norm(v-vm1)<epsilon
%                 break;
%           end         
    end
%     disp(k)
    
%     k
%     x0 = v;
     xRest = sign(v); 
     x0 = C*xRest; 
     x = 1/sqrt(2*par.N)*(x0(1:par.N,1)+1i*x0(par.N+1:2*par.N,1)); %*Q  *Q*par.N  1/sqrt(2*par.N)*
    
    beta = real(x'*H'*s)/(norm(H*x,2)^2+par.U*N0); 
%     beta = norm(x0)/sqrt(par.N*bps);
end




