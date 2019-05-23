% Bit  Error Rate 
%% developed by Lei Chu, Fei Wen and Robert Qiu
%% We develop this code based on https://github.com/quantizedmassivemimo/1bit_precoding;

%%% main simulation file for our published paper: Efficient Nonlinear Precoding for Massive
%%%  MU-MIMO Downlink Systems with 1-Bit DACs; 

%%%  notes:
%%%  1)   you should establish CVX which is avaliable at http://cvxr.com/cvx/
%%%  2)   if you use our codes, please cite our paper
%%%  3)   if you get any question for our codes, please contact    leochu@sjtu.edu.cn
%%%  

%%

close all; clear all; clc

parD.b = 1; % 1 means 1 bit
parD.U = 4; % number of UEs
parD.N =  16; % number of BS antennas
parD.trials = 1e3; % number of Monte-Carlo trials (transmissions)
parD.rHe = 0; % relative channel estimate error
parD.SNRdB_list = -10:2:14; % list of SNR [dB] values to be simulated 
parD.mod = 'QPSK'; % modulation type: 'QPSK','16QAM','64QAM'
parD.precoder =  {'MRT', 'WF', 'ZF','SQUID','Proposed','SDR1', 'ZFi'};

%% method_1

mod_methods = {'BPSK', 'QPSK', '8PSK', '16QAM','32QAM', '64QAM'};
mod = 'QPSK';

mod_order = find(ismember(mod_methods,mod));

if mod_order == 1
      mod_ind = 2^(mod_order-1);   
      sita = 0:pi/mod_ind:2*pi - pi/mod_ind;
      in_phase = cos(sita);
      quadrature = sin(sita);
      symbols = (in_phase + quadrature*1i)';
elseif mod_order == 2 || mod_order == 3 
      mod_ind = 2^(mod_order-1);   
      sita = 0:pi/mod_ind:2*pi - pi/mod_ind;
      in_phase = cos(sita + pi/4);
      quadrature = sin(sita + pi/4);
      symbols = (in_phase + quadrature*1i)';
      symbols = symbols/sqrt(mean(abs(symbols).^2));
else
      mod_ind = sqrt(2^(mod_order));   
      sita = 0:pi/mod_ind:2*pi - pi/mod_ind;
      in_phase = repmat(linspace(-1,1,mod_ind),mod_ind,1);
      quadrature = repmat(linspace(-1,1,mod_ind)',1,mod_ind);
      symbols = in_phase(:) + quadrature(:)*1i;   
      symbols = symbols/sqrt(mean(abs(symbols).^2));
end
parD.symbols = symbols';


%%

parD.card = length(parD.symbols);
parD.bps = log2(parD.card); 
parD.bits = de2bi(0:parD.card-1,parD.bps,'left-msb'); 



BER  = zeros(length(parD.b),length(parD.precoder),length(parD.SNRdB_list));
vrz = [];
for t=1:parD.trials
    t
    b = randi([0 1],parD.U,parD.bps);
   
    idx = bi2de(b,'left-msb')+1;
    s = parD.symbols(idx).';
%     s = parD.symbols(idx);
%     s = s';    
    
    n = sqrt(0.5)*(randn(parD.U,1)+1i*randn(parD.U,1));
    H = sqrt(0.5)*(randn(parD.U,parD.N)+1i*randn(parD.U,parD.N));
    H1 = sqrt(1 - parD.rHe)*H + ...
        sqrt(parD.rHe)*(randn(parD.U,parD.N)+1i*randn(parD.U,parD.N)); 
    for B = 1:length(parD.b)
            
             parD.quantizer = @(x) uqz(x, 1)/sqrt(parD.N); 
             parD.bg = sqrt(pi/2);   % bg decomposition factor
                         
            for pp=1:length(parD.precoder) 
                
                for k=1:length(parD.SNRdB_list)

                    N0 = 10.^(-parD.SNRdB_list(k)/10);
                    switch (parD.precoder{pp})
                        case 'MRTi'
                            [x, beta] = MRT(s,H1);
                        case 'MRT' 
                            [z, beta] = MRT(s,H1);
                            x = parD.quantizer(z); beta = beta/parD.bg;
                        case 'ZFi'
                            [x, beta] = ZF(s,H1); 
                        case 'ZF' 
                            [z, beta] = ZF(s, H1);
                            x = parD.quantizer(z); beta = beta/parD.bg;
                        case 'WFi' 
                            [x, beta] = WF(s,H1,N0);
                        case 'WF'  
                            [z, beta] = WF(s,H1,N0);
                            x = parD.quantizer(z); beta = beta/parD.bg;       
                        case 'SQUID' 
                            parD.b = 1; [x, beta,xRest] = SQUID(parD,s,H1,N0);  % only support for 1bit
                        case 'SDR1'    
                            parD.b = 1; 
                            [x, beta] = SDR_A1bit(parD,s,H1,N0);  % only support for 1bit
                        case 'Proposed' 
                            parD.b = 1;     
                            [x, beta, vr] = ADMM_Leo(parD,s,H1,N0);  vrz = [vrz vr];
                    end

                    Hx = H*x; 
                    y = Hx + sqrt(N0)*n;

                    shat = beta*y; 

                    [~,idxhat] = min(abs(shat*ones(1,length(parD.symbols))...
                            -ones(parD.U,1)*parD.symbols).^2,[],2); 
                    bhat = parD.bits(idxhat,:);

                    err = (idx~=idxhat); 
                    BER(B,pp,k) = BER(B,pp,k) + sum(sum(b~=bhat))/(parD.U*parD.bps);                   
                end         
            end     
    end
end

BER = squeeze(BER/parD.trials);


qq = -10:2:14;


semilogy(qq,BER(1,:),'-d',qq,BER(2,:),'-*',qq,BER(3,:) ,'-^',qq,BER(4,:) ,'-+',...
        qq,BER(5,:),'->',qq,BER(6,:),'-<',qq,BER(7,:) ,'-',...
        'LineWidth',2)

grid on


xlim([-10 14]);ylim([10^(-3) 10^(0)])
legend('MRT', 'WF', 'ZF',  'SQUID', 'Proposed', 'SDR', 'ZFi',3) 
xlabel('SNR (dB)')
ylabel('Uncoded BER')








