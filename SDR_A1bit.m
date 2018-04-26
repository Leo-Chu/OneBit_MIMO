function [x_sdr1, beta_sdr1] = SDR_A1bit(par,s,H,N0)

    HR0 = [ real(H) -imag(H) ; imag(H) real(H) ];
    sR = [real(s) ; imag(s)];
    
    if par.b == 1
        C = eye(2*par.N);  bps = 2; 
        HR = HR0*C; num = bps/2; N0 = N0*num^2; % 这里没有 虚拟向量扩充，所以信噪比没有变化！！
        cvx_begin quiet
            TR = [HR'*HR+par.U*N0*eye(2*par.N), -HR'*sR; -sR'*HR, norm(sR,2)^2];
            variable BR(2*par.N+1,2*par.N+1) symmetric
            minimize trace(TR*BR)
            subject to
            for j = 2:2*par.N 
                BR(1,1) == BR(j,j);
            end
            BR(2*par.N+1,2*par.N+1) == 1;
            BR == semidefinite(2*par.N+1);
        cvx_end
    elseif par.b == 2
        Q = 20*log10(5);    
        C = [2*eye(2*par.N) eye(2*par.N)]; bps = 4; num = bps/2; N0 = N0*Q;
        HR = HR0*C;% size(HR)
        cvx_begin quiet
            TR = [HR'*HR+par.U*N0*eye(2*num*par.N), -HR'*sR; -sR'*HR, norm(sR,2)^2];
            variable BR(2*num*par.N+1,2*num*par.N+1) symmetric
            minimize trace(TR*BR)
            subject to
            for j = 2:2*num*par.N 
                BR(1,1) == BR(j,j);
            end
            BR(2*num*par.N+1,2*num*par.N+1) == 1;
            BR == semidefinite(2*num*par.N+1);
        cvx_end
    elseif par.b == 3
        Q = 20*log10(21);
        C = [4*eye(2*par.N) 2*eye(2*par.N) eye(2*par.N)]; bps = 6; num = bps/2; N0 = N0*Q;
        HR = HR0*C;% size(HR)
        cvx_begin quiet
            TR = [HR'*HR+par.U*N0*eye(2*num*par.N), -HR'*sR; -sR'*HR, norm(sR,2)^2];
            variable BR(2*num*par.N+1,2*num*par.N+1) symmetric
            minimize trace(TR*BR)
            subject to
            for j = 2:2*num*par.N 
                BR(1,1) == BR(j,j);
            end
            BR(2*num*par.N+1,2*num*par.N+1) == 1;
            BR == semidefinite(2*num*par.N+1);
        cvx_end
    else
        disp('Not supportted at current version!!!');
    end

    MAX = 1e2;
    
    % eigenvalue decomposition
    [V, D] = eig(BR);
    
    % -- Rank-one approximation
    
    % find maximum eigenvalue 
    [~, idxmax] = max(diag(D)); 
    
    % quantize to feasible solution
    xR_sdr11 = par.quantizer(sqrt(D(idxmax,idxmax))*V(:,idxmax));

%     xR_sdr11 = sign(sqrt(D(idxmax,idxmax))*V(:,idxmax));
    xR_sdr1 = C*xR_sdr11(1:2*num*par.N);
    
    x_sdr1 = sign(xR_sdr1(end)) * (xR_sdr1(1:par.N,1)+1i*xR_sdr1(par.N+1:2*par.N,1));
    %
%     x_sdr1 =  C*x_sdr11; % duo bits
        
    % compute precoding factor
    beta_sdr1 = real(x_sdr1'*H'*s)/(norm(H*x_sdr1,2)^2+par.U*N0);
    
end

