% convergence analysis

close all; clear all; clc

parD.b = 1; % 1 means 1 bit
parD.U = 4; % number of UEs
parD.N =  16; % number of BS antennas
parD.trials = 1e0; % number of Monte-Carlo trials (transmissions)
parD.rHe = 0; % relative channel estimate error
parD.SNRdB_list = -5:5:15; % list of SNR [dB] values to be simulated 
parD.mod = 'QPSK'; % modulation type: 'QPSK','16QAM','64QAM'
parD.precoder =  {'ADMM_Leo1'};


%% method_1

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

parD.card = length(parD.symbols);
parD.bps = log2(parD.card); 
parD.bits = de2bi(0:parD.card-1,parD.bps,'left-msb'); 

load('codebook_downlink.mat') 


BER  = zeros(length(parD.b),length(parD.precoder),length(parD.SNRdB_list));
 vrz2 = [];    wrz2 = [];  
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
        sqrt(parD.rHe/2)*(randn(parD.U,parD.N)+1i*randn(parD.U,parD.N)); 

    for B = 1:length(parD.b)
            
            parD.B = 2^B; 
            parD.Q = log2(parD.B); 
            parD.lsb = lsb_list(parD.B-1)/sqrt(2*parD.N); 
%             [~, parD.labels, parD.thresholds, delta] = uqz(1, 2^parD.b ); 
             parD.quantizer = @(x) uqz(x, 1)/sqrt(parD.N);    
%              parD.quantizer = @(x) uqz(x, parD.lsb, 1)/sqrt(parD.N); 

%              parD.labels = parD.labels; 

             parD.bussgang = parD.lsb*sqrt(parD.N/pi)...
                *sum(exp(-parD.N*parD.lsb^2*((1:parD.B-1)-parD.B/2).^2)); 
            for pp=1:length(parD.precoder) 
                 vrz1 = [];  wrz1 = []; 
                % SNR loop
                for k=1:length(parD.SNRdB_list)

                    N0 = 10.^(-parD.SNRdB_list(k)/10);
                    switch (parD.precoder{pp})

                        case 'ADMM_Leo1' 
                            parD.b = 1;     
                            [x, beta,wr, vr] = ADMM_Leo(parD,s,H1,N0);  
                            
%                         case 'ADMM_Leo1x' 
%                             parD.b = 1;     
%                             [x, beta,wr, vr] = ADMM_Leo1x(parD,s,H1,N0);  
%                             vrz2 = [vrz2 vr]; wrz2 = [wrz2 wr]; 
                    end   
                    vrz1 = [vrz1; vr]; wrz1 = [wrz1; wr]; 
                end 
                    vrz2 = [vrz2; vrz1]; wrz2 = [wrz2; wrz1]; 
            end     
    end
end

qq = 1:size(vrz1,2);

wrz2(find(wrz2==0))=7.850462293418876e-17;

semilogy(qq,wrz2(1,:),'-',qq,wrz2(2,:),'-',qq,wrz2(3,:) ,'-',qq,vrz2(4,:) ,'-', qq,vrz2(5,:),'-',...
        'LineWidth',2) % ,'-^',qq,vrz2(4,:) ,'-+',... qq,vrz2(5,:),'-'

grid on

legend('SNR = -5 dB', 'SNR = 0 dB', ...
        'SNR = 5 dB','SNR = 10 dB','SNR = 15 dB') 
xlabel('Iterations')
ylim([10^(-20) 10^(0)]) 
ylabel({'$$\Delta {\tilde{\bf{v}}^{k}}$$'},'Interpreter','latex')

%%

% qq = 1:size(vrz1,2);
% 
% % wrz2(find(wrz2==0))=7.850462293418876e-17;
% 
% semilogy(qq,vrz2(1,:),'-',qq,vrz2(2,:),'-',qq,vrz2(3,:) ,'-',qq,vrz2(4,:) ,'-',qq,vrz2(5,:),'-',...
%         'LineWidth',2) % ,'-^',qq,vrz2(4,:) ,'-+',... qq,vrz2(5,:),'-'
% 
% grid on
% 
% % set(gca,'FontSize',14);  
% legend('SNR = -5 dB', 'SNR = 0 dB', ...
%         'SNR = 5 dB','SNR = 10 dB','SNR = 15 dB') 
% xlabel('Iterations')
% ylabel({'$$\Delta {{\bf{u}}^{k}}$$'},'Interpreter','latex')


