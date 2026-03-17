function Iuv = paperDiffuse(DIMS,L,h,g,varargin)
DIMS = [2048,2048];

h2 = 10; % Distance of target from card in mm 
         % (not equivalent to working distance)

u = linspace(-1,1,DIMS(1));
v = linspace(-1,1,DIMS(2));
[U,V] = ndgrid(u,v);
uu = real(sqrt(U.^2+V.^2));
COS_TH = sqrt(1-(U.^2+V.^2));
COS_TH(abs(imag(COS_TH))>0) = 0;
SEC_TH = 1./COS_TH.*(COS_TH>0);
SIN_TH = sqrt(U.^2+V.^2);
vv = U.*SEC_TH.*SIN_TH;
uu(abs(imag(uu))>0) = 0;
vv(abs(imag(vv))>0) = 0;

Iuv = h.*COS_TH.^3.*(h.*cos(g)+L.*sin(g)-h2.*vv.*sin(g))./(h2.^2.*((h.^2+L.^2)+2.*h2.*vv.*(h2.*vv-L))).^2;
Iuv(isnan(Iuv)) = 0;
end
