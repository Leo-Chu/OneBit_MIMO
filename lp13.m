% differ modulations mode

close all; clear all; clc

parD.b = 1; % 1 means 1 bit
parD.U = 10; % number of UEs
parD.N =  128; % number of BS antennas
parD.trials = 2e3; % number of Monte-Carlo trials (transmissions)
parD.rHe = 0; % relative channel estimate error
parD.SNRdB_list = -10:2:14; % list of SNR [dB] values to be simulated 
mod = {'QPSK','16QAM','64QAM'}; % modulation type: 'QPSK','16QAM','64QAM'
parD.precoder =  {'ADMM_Leo1','ZFi'};

%%

BER  = zeros(length(mod),length(parD.precoder),length(parD.SNRdB_list));
vrz = [];
for t=1:parD.trials
    t
    for m = 1:length(mod)
            parD.mod = mod{m};
            
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
            parD.card = length(parD.symbols);
            parD.bps = log2(parD.card); 
            parD.bits = de2bi(0:parD.card-1,parD.bps,'left-msb'); 
            
            b = randi([0 1],parD.U,parD.bps);
   
            idx = bi2de(b,'left-msb')+1;
            s = parD.symbols(idx).';

            n = sqrt(0.5)*(randn(parD.U,1)+1i*randn(parD.U,1));
            H = sqrt(0.5)*(randn(parD.U,parD.N)+1i*randn(parD.U,parD.N));
            H1 = sqrt(1 - parD.rHe)*H + ...
                sqrt(parD.rHe)*(randn(parD.U,parD.N)+1i*randn(parD.U,parD.N)); 
            
            for pp=1:length(parD.precoder) 

                % SNR loop
                for k=1:length(parD.SNRdB_list)

                    N0 = 10.^(-parD.SNRdB_list(k)/10);
                    switch (parD.precoder{pp})
                        case 'ZFi'
                            [x, beta] = ZF(s,H1); 
                        case 'ADMM1' 
                            parD.b = 1; [x, beta,xRest] = ADMM_Mbits(parD,s,H1,N0);  
                        case 'ADMM_Leo1' 
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
                    BER(m,pp,k) = BER(m,pp,k) + sum(sum(b~=bhat))/(parD.U*parD.bps);                   
                end         
            end     
    end
end

BER = squeeze(BER/parD.trials);


qq = -10:2:14;

% semilogy(qq,BER(1,:),'-k.',qq,BER(2,:),'-g*',qq,BER(3,:),'-b^',qq,BER(4,:),'-c+',...
%         qq,BER(5,:),'-m^','LineWidth',2)  %qq,BER(6,:),'-y+',
% legend('MRT', 'WF', 'ZF',  'SQUID', 'ZFi',3) 


p=semilogy(qq,squeeze(BER(1,1,:)),'-bd',qq,squeeze(BER(2,1,:)),'-bo',...
        qq,squeeze(BER(3,1,:)) ,'-b+',qq,squeeze(BER(1,2,:)) ,'-kd',...
        qq,squeeze(BER(2,2,:)),'-ko',qq,squeeze(BER(3,2,:)),'-k+',...
        'LineWidth',1.7);
% , qq,BER(8,:) ,'-k',
grid on

% set(gca,'FontSize',14);
xlim([-10 14]); ylim([10^(-4) 10^(0)])
% set(gca,'xaxislocation','top');  set(gca,'yaxislocation','right');
% legend('PP-QPSK', 'PP-16QAM',...
%          'PP-64QAM', ...
%        'ZFi-QPSK', 'ZFi-16QAM','ZFi-64QAM',4) 

% legend(p(1:3),'Proposed-QPSK','Proposed-16QAM','Proposed-64QAM');
% ah=axes('position',get(gca,'position'),...
%             'visible','off'); set(gca,'FontSize',14);
% legend(ah,p(4:6),'ZFi-QPSK','ZFi-16QAM','ZFi-64QAM','location','west');

legend('Proposed-QPSK','Proposed-16QAM','Proposed-64QAM',...
        'ZFi-QPSK','ZFi-16QAM','ZFi-64QAM')
xlabel('SNR (dB)')
ylabel('Uncoded BER')