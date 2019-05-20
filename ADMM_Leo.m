%% developed by Lei Chu and Fei Wen

function [x, beta, zr,vr] = ADMM_Leo(par,s,H,N0)

    HR0 = [ real(H) -imag(H) ; imag(H) real(H) ];
    sR = [real(s) ; imag(s)];
    
    if par.b == 1
        C = eye(2*par.N);  bps = 2; 
    elseif par.b == 2
        Q1 = sqrt(5);    
        C = [2*eye(2*par.N) eye(2*par.N)]/Q1; bps = 4; 
    elseif par.b == 3
        Q1 = sqrt(21);

        C = [4*eye(2*par.N) 2*eye(2*par.N) eye(2*par.N)]/Q1; bps = 6; 
    else
        disp('Not supportted at current version!!!');
    end
    HR = HR0*C;

    v = zeros(par.N*bps,1);
    z = zeros(par.N*bps,1);
    w = zeros(par.N*bps,1);
    
    iter = 5e2; c = par.U*N0; epsilon = 1e-6;
   vr=[]; zr = [];
   rho_0=1;
   rho = rho_0;
   rho_t = 2.1*norm(HR'*HR,2); 
   [uu,ss,~] = svd(HR'*HR);
   Hs = HR'*sR;

    for k = 1:iter
          vm1 = v;
          zm1 = z; 
          if rho<rho_t
                rho = rho_0*1.15^k;
          end
          gg = 2*diag(ss) + (2*c+rho);

          v = uu*((uu'*(2*Hs + rho*z + w))./gg);
          z = sign(v - w/rho)*norm(v - w/rho,1)/(par.N*bps);
          w = w - rho*(v - z);
          
          vr = [vr, norm(v-vm1)/norm(v)];
          zr = [zr, norm(z-zm1)/norm(z)];
          
%  line 50-52 should be deleted when peforming convergence analysis!!!!!!          
          if norm(v-vm1)<epsilon
                break;
          end         
    end

     xRest = sign(v); 
     x0 = C*xRest; 
     x = 1/sqrt(2*par.N)*(x0(1:par.N,1)+1i*x0(par.N+1:2*par.N,1)); 
    
    beta = real(x'*H'*s)/(norm(H*x,2)^2+par.U*N0); 

end




