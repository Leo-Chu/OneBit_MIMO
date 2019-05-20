

function [output, labels, tau, clip] = uqz(y, lsb, L)
        if(L==1)             
                output = (sign(real(y)) + 1i*sign(imag(y)))/sqrt(2);
                labels = 1;   tau = 1;   clip = 1;          
        else
                clip = lsb*L/2;

                % clip signal
                if isreal(y)
                yc = max(min(y,clip-lsb/1e5),-(clip-lsb/1e5));
                else
                yc = max(min(real(y),clip-lsb/1e5),-(clip-lsb/1e5)) + 1i*max(min(imag(y),clip-lsb/1e5),-(clip-lsb/1e5));
                end

                % quantizer
                if mod(L,2) == 0
                Q = @(x) lsb*floor(x/lsb) + lsb/2; % midrise quantizer (without clipping)
                else
                Q = @(x) lsb*floor(x/lsb + 1/2); % midtread quantizer (without clipping)
                end

                % quantize signal
                if isreal(y)
                output = Q(yc);
                else
                output = Q(real(yc)) + 1i*Q(imag(yc));
                end


                % uniform quantization labels
                labels = lsb *((0:L-1) - (L-1)/2);

                % uniform quantization thresholds
                tau = [-10^100, bsxfun(@minus, labels(:,2:end), lsb/2), 10^100];	
        end
        
end








