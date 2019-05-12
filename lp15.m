% differ computational time

close all; clear all; clc

parD.b = 1; % 1 means 1 bit
parD.U = 4; % number of UEs
NT =  [16, 32, 64, 128, 256]; % number of BS antennas
parD.trials = 1e3; % number of Monte-Carlo trials (transmissions)
parD.rHe  = 0;% relative channel estimate error
parD.SNRdB_list = 0; % list of SNR [dB] values to be simulated 
parD.mod = 'QPSK'; % modulation type: 'QPSK','16QAM','64QAM'
parD.precoder =  {'MRT', 'WF', 'ZF','SQUID','ADMM','SDR1'};


switch (parD.mod)
    case 'QPSK'
        parD.symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ];
    case '16QAM'
        parD.symbols = [...
            -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM'
        parD.symbols = [...
            -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
end

parD.E =  sqrt(mean(abs(parD.symbols).^2));
parD.symbols = parD.symbols/parD.E;
% scatterplot(symbol_bolck)

%%
N0 = 10.^(-parD.SNRdB_list/10);
parD.card = length(parD.symbols);
parD.bps = log2(parD.card); 
parD.bits = de2bi(0:parD.card-1,parD.bps,'left-msb'); 

load('codebook_downlink.mat') 


CT  = zeros(length(parD.precoder),length(NT));
vrz = [];
for t=1:parD.trials
    t
    b = randi([0 1],parD.U,parD.bps);
   
    idx = bi2de(b,'left-msb')+1;
    s = parD.symbols(idx).';
    
    n = sqrt(0.5)*(randn(parD.U,1)+1i*randn(parD.U,1));
    
    parD.B = 2^parD.b; 
    parD.Q = log2(parD.B); 

    for pp=1:length(parD.precoder) 

        for k=1:length(NT)
            parD.N = NT(k);    
            parD.lsb = lsb_list(parD.B-1)/sqrt(2*NT(k)); 
            parD.quantizer = @(x) uqz(x,1)/sqrt(NT(k)); 
            parD.bussgang = parD.lsb*sqrt(NT(k)/pi)...
                *sum(exp(-NT(k)*parD.lsb^2*((1:parD.B-1)-parD.B/2).^2)); 
            H = sqrt(0.5)*(randn(parD.U,NT(k))+1i*randn(parD.U,NT(k)));
            H1 = sqrt(1 - parD.rHe)*H + ...
                sqrt(parD.rHe/2)*(randn(parD.U,NT(k))+1i*randn(parD.U,NT(k))); 
            tt = 0;
            switch (parD.precoder{pp})

                case 'MRT' 
                    tic
                    [z, beta] = MRT(s,H1);
                    tt = toc;
                    x = parD.quantizer(z); beta = beta/parD.bussgang;
                case 'ZFi'
                    tic 
                    [x, beta] = ZF(s,H1);  tt = toc;
                case 'ZF' 
                    tic 
                    [z, beta] = ZF(s, H1); tt = toc;
                    x = parD.quantizer(z); beta = beta/parD.bussgang;
                case 'WF'  
                    tic
                    [z, beta] = WF(s,H1,N0); tt = toc;
                    x = parD.quantizer(z); beta = beta/parD.bussgang;       
                case 'SQUID'       
                    parD.b = 1; tic
                    [x, beta,xRest] = SQUID(parD,s,H1,N0);  tt = toc;
                case 'SDR'    
                    parD.b = 1; tic 
                    [x, beta] = SDR_A1bit(parD,s,H1,N0); tt = toc;
                case 'ADMM' 
                    parD.b = 1;   tic  
                    [x, beta, ~,~] = ADMM_Leo(parD,s,H1,N0); tt = toc;% vrz = [vrz vr];
            end

            Hx = H*x; 
            y = Hx + sqrt(N0)*n;
            CT(pp,k) = CT(pp,k) + tt;    
        end           
    end
end

qq = [16, 32, 64, 128, 256];
% CT = CT/parD.trials; % computaional time

semilogy(qq,CT(1,:),'-d',qq,CT(2,:),'-*',qq,CT(3,:) ,'-^',qq,CT(4,:) ,'-+',...
        qq,CT(5,:),'->',qq,CT(6,:),'-<',qq,CT(7,:) ,'-',...
        'LineWidth',2);

grid on


xlim([16 256])

legend('MRT', 'WF', 'ZF',  'SQUID', 'Proposed', 'SDR', 'ZFi',3) 

xlabel('Number of transmit antennas')
ylabel('Computaional time')
