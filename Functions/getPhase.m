% Paloma Casteleiro - 09.02.18

function [Phase_im, IPhase] = getPhase(im_in, HF,type)
% Function gets phase image from DPC image using TIAN15 method

dim = size(im_in{1});
Inum = zeros(size(im_in{1}));
Iden = zeros(size(im_in{1}));

for col=1:2
    % Filter Image
    I = im_in{col};
    
    % Deconvolve Image from TF to get Quantitative Phase
    Hf = HF{col}.*(-1)^(col+1);
    If = fftshift(fft2(I));   
    if strcmp(string(type),"Paper")
        a = 2e-2;
    elseif strcmp(string(type),"Brain")
        a = 5e-3;
    end
    
    Ih1 = ifft2(ifftshift((-1i*Hf.*If)./(abs(Hf).^2+a)));
    IPhase{col} = real(Ih1);
    
    Inum = Inum+(-1i*Hf.*If);
    Iden = Iden+(abs(Hf).^2+a);
end

Phase_im = real(ifft2(ifftshift(Inum./Iden)));

end
