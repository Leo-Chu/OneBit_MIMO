

function [output] = uqz(y, L)
        if(L==1)             
                output = (sign(real(y)) + 1i*sign(imag(y)))/sqrt(2);       
        else
			print('not supported now!')
        end
        
end








