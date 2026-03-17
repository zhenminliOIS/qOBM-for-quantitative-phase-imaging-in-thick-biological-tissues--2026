function [filtR]=RadialFiltGen(dim1,dim2,filter_type,FreqRange,filt_pow)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [filtY]=FiltGen(dim1,dim2,filtdim,filter_type,center,width,twosided,filt_pow)
% 
% INPUTS: 
%        dim1: number of rows
%        dim2: number of columns
%        filter type:
%            1: low pass filter
%            2: butterworth filter 
%           -2: butterworth filter high pass filter
%       FreqRange: cutoff of filter in pixles
%       filt_pow: specify power of filter-only for butterworth            
%
% OUTPUTS: 
%       Filter matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[X, Y] = meshgrid((1:dim2)- dim2/2,(1:dim1)-dim1/2);
R = (X.^2 + Y.^2).^(0.5);

if filter_type == 1 % Low pass filter hard cut off
    filtR = 1.*R<FreqRange;
        
elseif abs(filter_type) == 2 % Butterworth filter; negative if high pass
   
        filtR = (1+((R)/(FreqRange/2)).^(2*filt_pow)).^(-1/2);
        if filter_type<0; filtR = 1-filtR; end
    
end

