function [TF,B,varargout] = get3DTFnpbb(X,Y,Z,S,P,lds,varargin)
% get3DTFnpbb Patrick Ledwig 12/8/21
%   get the rapid broadband optical transfer function for transmission
%   optical system
% INPUTS:
%   X,Y,Z - Meshgrid variables with spatial coordinates. If
%       options.KSPACE==1, then these coordinates are given by U,V,W in
%       units of m^-1/nu_max
%   S - Source angular distribution in u,v from 0 to 1 where
%       u=sin(th)cos(phi), v=sin(th)sin(phi), where th and phi are the
%       axial and azimuthal angle of incidence. if this is an NxM matrix,
%       then the angular distribution is taken to be uniform over the
%       source temporal spectrum. if this is an NxMxL matrix, where L is
%       any integer >=2 then the source is interpolated over the lambda
%       range. The values are relative. If this is a single value, it will
%       be taken to be the source aperture NA.
%   G - Source relative temporal spectrum. Lx1 array. If L!=numel(lda),
%       then it is interpolated over the provided lambdas. For an arbitrary
%       distribution, make S NuxNvxNw, with G=1
%   lda - Source temporal spectrum abscissa (wavelength or c-normalized
%       frequency). If contains one value, single wavelength code will run.
%   P - microscope pupil function. If this is a single number, it will be
%       taken as the NA. If it is NxMx1, then phase will be adjusted
%       according to wavelength, with the phase provided applying to nu_max
%   options - text flags for the operation of the program
%       KSPACE - boolean flag determining the units of the input coords
%       NU - boolean flag, if set, lda is given in c-normalized freq
%       GPU - boolean flag allowing CUDA computation
%       Nr - number of points to sample in radius on Ewald plane [15]
%       Nph - number of points to sample in angle on Ewald plane [5]
%       verbose - flag to output to console
%       plot - flag to plot results

% EXAMPLE VARIABLES
% P=.8;
% Nx=127; Ny=127; Nz=63; Nl=10;
% dx=630e-3/4;
% dy=630e-3/4;
% dz=630e-3/2;
% x=((0:(Nx-1))-floor(Nx/2))*dx;
% y=((0:(Ny-1))-floor(Ny/2))*dy;
% z=((0:(Nz-1))-floor(Nz/2))*dz;
% [X,Y,Z]=ndgrid(x,y,z);
% uMax=1/dx;
% vMax=1/dy;
% u=linspace(-1,1,Nx); v=linspace(-1,1,Ny);
% lmax=.800;
% lmin=.44;
% Nr0=100;
% lds=linspace(lmin,lmax,Nr0);
% [U,V]=ndgrid(u,v);
% G=ones(Nl,1);
% S=double(sqrt(U.^2+V.^2)<=.65);
% S=S.*(U>0);
% % S=.4;
% varargin={'plot','verbose','GPU','Nr',20,'Nph',15};

% q1i=linspace(-2,2,Nx);
% q2i=linspace(-2,2,Ny);
% q3i=linspace(-1,1,Nz);
% if Nx==1 q1i=-0; end
% if Ny==1 q2i=-0; end
% if Nz==1 q3i=-0; end
% [q1,q2,q3]=ndgrid(q1i,q2i,q3i);

% Set flags
if length(varargin)>0
    options=varargin{1};
end
if length(varargin)>1
    updims=varargin{2};
else
    updims=[255,255,63];
end
if any(strcmpi(options,'GPU'))
    GPUFLAG=1;
else
    GPUFLAG=0;
end
if any(strcmpi(options,'Nbatch'))
    Nbatch=Nbatch;
else
    if GPUFLAG
        gpuBlock=gpuDevice().MaxThreadBlockSize;
        Nbatch=10*prod(ceil(size(X)./gpuBlock));
    else
        Nbatch=100;
    end
end
if any(strcmpi(options,'verbose'))
    VERBOSE=1;
else
    VERBOSE=0;
end
if any(strcmpi(options,'GPU'))
    GPU=true;
else
    GPU=false;
end
if any(strcmpi(options,'G'))
    Nr=options{strcmpi(options,'G')+1};
else
    G=1;
end
if any(strcmpi(options,'Nr'))
    Nr=options{strcmpi(options,'Nr')+1};
else
    Nr=15;
end
if any(strcmpi(options,'Nph'))
    Nph=options{strcmpi(options,'Nph')+1};
else
    Nr=5;
end
if any(strcmpi(options,'NU'))
    NU=true;
else
    NU=false;
end
if any(strcmpi(options,'Born'))
    BORN=true;
else
    BORN=false;
end

if numel(P)==1
    NAp=P;
    PFLAG=false;
else
    NAp=1;
    PFLAG=true;
end
if numel(S)==1
    NAs=S;
else
    NAs=NAp;
end
if isreal(P)
    PCX=false;
else
    PCX=true;
end
B=1;

myFastTime=tic;
tic;
lambda=(min(lds)+max(lds))/2;
%%%%
[Nx,Ny,Nz]=size(X);
%Define image space
x=sort(unique(X));
y=sort(unique(Y));
z=sort(unique(Z));
dx=x(2)-x(1);
dy=y(2)-y(1);
dz=z(2)-z(1);

uMax=1/dx;
vMax=1/dy;
wMax=1/dz;

%Define the down-sample space for the output
ui=((0:(Nx-1))-floor(Nx/2))*uMax/Nx*lambda;
vi=((0:(Ny-1))-floor(Ny/2))*vMax/Ny*lambda;
wi=((0:(Nz-1))-floor(Nz/2))*wMax/Nz*lambda;
[Ui,Vi,Wi]=ndgrid(ui,vi,wi);
[Ui2,Vi2]=ndgrid(ui,vi);
dui=ui(2)-ui(1);
dvi=vi(2)-vi(1);
dwi=wi(2)-wi(1);

if any(strcmpi(options,'upsample'))
    Nx0=updims(1); Ny0=updims(2); Nz0=updims(3);
    x0=((0:(Nx0-1))-floor(Nx0/2))*dx;
    y0=((0:(Ny0-1))-floor(Ny0/2))*dy;
    z0=((0:(Nz0-1))-floor(Nz0/2))*dz;

    %Define un-normalized 2d angle space for upsampling
    u=((0:(Nx0-1))-floor(Nx0/2))*uMax/Nx0;
    v=((0:(Ny0-1))-floor(Ny0/2))*vMax/Ny0;
    w=((0:(Nz0-1))-floor(Nz0/2))*wMax/Nz0;

    Nx=Nx0; Ny=Ny0; Nz=Nz0;
else
    u=((0:(Nx-1))-floor(Nx/2))*uMax/Nx;
    v=((0:(Ny-1))-floor(Ny/2))*vMax/Ny;
    w=((0:(Nz-1))-floor(Nz/2))*wMax/Nz;
end



%Define normalized 2D angular space from image space
ul=u*lambda; vl=v*lambda;
dul=ul(2)-ul(1); dvl=vl(2)-vl(1);
[Ul,Vl]=ndgrid(ul,vl);
%Define 3rd dimension for 3D angle space
wl=w*lambda; dwl=wl(2)-wl(1);
[Ul3,Vl3,Wl3]=ndgrid(ul,vl,wl);
%Assign a pupil space
if PFLAG
    R=sqrt(Ul.^2+Vl.^2);
    if PCX
        Pang=angle(P);
        P=abs(P);
        Fpang=griddedInterpolant(Ui2,Vi2,Pang,'linear','none');
        Pi_ang=Fpang(Ul,Vl); Pi_ang(isnan(Pi_ang))=0;
    end
    Fp=griddedInterpolant(Ui2,Vi2,P,'linear','none');
    Pi=Fp(Ul,Vl); Pi(isnan(Pi))=0;
    indPi=find(Pi~=0);
    NAp=max(R(indPi(:)));

end
%%%%

q1=Ul3; q2=Vl3; q3=Wl3;
qx=sqrt(q1.^2+q2.^2);
qx=abs(qx);
alph=tan(asin(NAs));
bet=tan(asin(NAp));
q=sqrt(qx.^2+q3.^2);
ca=1/sqrt(1+alph.^2);
cb=1/sqrt(1+bet.^2);
% GPU=true;
tic;

%Establish Source
if numel(S)==1
    S0=(sqrt(U0.^2+V0.^2)<=S);
    Nu0=128;Nv0=128;Nr0=length(G);
elseif numel(size(S))==2
    [Nu0,Nv0]=size(S);
    Nr0=length(G);
    S0=S;
elseif numel(size(S))==3
    [Nu0,Nv0,Nr0]=size(S);
    S0=S;
end
u0=linspace(-1,1,Nu0);
v0=linspace(-1,1,Nv0);
du0=u0(2)-u0(1); dv0=v0(2)-v0(1);
[U0,V0]=ndgrid(u0,v0);
R0=sqrt((U0-.5).^2+V0.^2);
R1=sqrt((U0+.5).^2+V0.^2);
COS0=real(sqrt(1-U0.^2-V0.^2));
S0=S0.*sum(abs(S0(:))>0)./sum(abs(S0(:)));


if numel(lds)==1 %single value of G defined
    dlds=sqrt(dul.^2+dvl.^2+dwl.^2)/2;
    lmin=lds-dlds;
    lmax=lds+dlds;
    vmin=lmin./lmax;
    vmax=1;
    nu0=linspace(vmin,vmax,3);
    dnu0=vmax-vmin;
    G=[.25;.5;.25];
else %multiple values of lambda defined
    lmin=min(lds);
    lmax=max(lds);
    vmin=min(lds)./max(lds);
    vmax=1;
    nu0=linspace(vmin,vmax,Nr0);
    dnu0=nu0(2)-nu0(1);
end

G=reshape(nu0,1,1,length(nu0));
[Uf,Vf,Rf]=ndgrid(u0,v0,nu0);
cosf=sqrt(1-Uf.^2-Vf.^2);
Nr0=length(nu0);
Sf0=repmat(S0,1,1,Nr0).*repmat(G,Nu0,Nv0,1);
Sf0sum=sum(abs(Sf0(:)).*du0.*dv0.*dnu0);
Sf=repmat(S,1,1,Nr0).*repmat(G,Nu0,Nv0,1)./Sf0sum;
Fsrc=griddedInterpolant(Uf,Vf,Rf,Sf,'linear','none');

if PFLAG %For arbitrary input Pupil function
    P0=Fp(U0,V0); P0(isnan(P0))=0;
    Pf0=repmat(P0,1,1,Nr0);
    Pf0sum=sum(abs(Pf0(:)).*du0.*dv0.*dnu0);
    Pf=repmat(P0,1,1,Nr0)./Pf0sum;
    if PCX
        P0_ang=Fpang(U0,V0); P0a_ng(isnan(P0_ang))=0;
        Pf0_ang=repmat(P0_ang,1,1,Nr0);
        Pf0_angsum=sum(abs(Pf0_ang(:)).*du0.*dv0.*dnu0);
        Pf_ang=repmat(P0_ang,1,1,Nr0);
        Fppl_ang=griddedInterpolant(Uf,Vf,Rf,Pf_ang,'linear','none');  
    end
     Fppl=griddedInterpolant(Uf,Vf,Rf,Pf,'linear','none');
end

Nq=numel(q);
qInd=reshape(1:numel(q),Nx,Ny,Nz);
Tab=squeeze(zeros(size(q)));
Tph=squeeze(zeros(size(q)));
qMask=zeros(1,1,Nq);

Nr=15;
Nph=5;
Nq=numel(q);
qInd=reshape(1:numel(q),Nx,Ny,Nz);
Tab=squeeze(zeros(size(q)));
Tph=squeeze(zeros(size(q)));
qMask=zeros(1,1,Nq);

totalTime=tic;
% tc=[68,64,56]; %test point for debugging
Nrk=2;
Nphk=2;
for kk=1:2
    sgn=(-1)^(kk+1);
    q=reshape(q,1,1,Nq);
    q3=reshape(q3,1,1,Nq);
    qx=reshape(qx,1,1,Nq);
    rhmin=real(sqrt(complex(vmin.^2-q.^2/4)));
    rhmax=real(sqrt(complex(vmax.^2-q.^2/4)));
    rh=[1:Nrk]'./Nrk.*(rhmax-rhmin)+rhmin;
    cosphs=-real(q./(qx.*rh).*(-sgn*q3/2+ca.*sqrt(rh.^2+q.^2./4)));
    cosphp=-real(q./(qx.*rh).*(sgn*q3/2+cb.*sqrt(rh.^2+q.^2./4)));
    cosmax= min((cosphs),(cosphp));
    phix=mod(real(acos(complex(cosmax))),pi);
    phixmin= phix;
    phixmax= mod(2*pi-phix,2*pi);
    dph=sum((phixmax-phixmin),1)./(Nphk*Nrk);
    qMask=qMask|(dph>0);
end
mInd=find(qMask);
% myInd=sub2ind([Nx,Ny,Nz],tc(1),tc(2),tc(3)); %test point for debugging
% myMInd=find(mInd==myInd);
Nm=numel(mInd);
q=q(mInd);
qx=qx(mInd);
q1=q1(mInd);
q2=q2(mInd);
q3=q3(mInd);
tab=Tab(mInd);
tph=Tph(mInd);
qmInd=1:numel(q);
if GPU
    q=gpuArray(q);
    q1=gpuArray(q1);
    q2=gpuArray(q2);
    q3=gpuArray(q3);
    qx=gpuArray(qx);
    tab=gpuArray(tab);
    tph=gpuArray(tph);
end
if BORN&&PFLAG&&PCX
    tab=complex(tab);
    tph=complex(tph);
end

for bb=1:Nbatch
    batchTime=tic;
    bInd=qmInd((bb-1)*floor(Nm./Nbatch)+1:(bb)*floor(Nm./Nbatch));
    Nqb=numel(bInd);
    qb=reshape(q(bInd),1,1,Nqb);
    q1b=reshape(q1(bInd),1,1,Nqb);
    q2b=reshape(q2(bInd),1,1,Nqb);
    q3b=reshape(q3(bInd),1,1,Nqb);
    qxb=reshape(qx(bInd),1,1,Nqb);

    rhmin=real(sqrt(complex(vmin.^2-qb.^2/4)));
    rhmax=real(sqrt(complex(vmax.^2-qb.^2/4)));
    rh=[1:Nr]'./Nr.*(rhmax-rhmin)+rhmin;
    for kk=1:2 %for each slice
        sgn=(-1)^(kk+1);
        cosmax=real(-qb./(2*qxb.*rh).*(sqrt(rh.^2+qb.^2/4).*(ca+cb)+abs(sqrt(rh.^2+qb.^2/4).*(ca-cb)-sgn.*q3b)));
        phix=real(acos(complex(cosmax)));
        phixmin= phix;
        phixmax= mod(2*pi-phix,2*pi);
        phi=[0.5:Nph-.5]./Nph.*(phixmax-phixmin)+phixmin;
        s1= rh.*cos(phi);
        s2= rh.*sin(phi);
        u1= sgn.*q1b/2-q2b.*s2./qxb+q1b.*q3b.*s1./(qxb.*qb);
        u2= sgn.*q2b/2+q1b.*s2./qxb+q2b.*q3b.*s1./(qxb.*qb);
        u3= sgn.*q3b/2-qxb.*s1./qb;
        rf= sqrt(u1.^2+u2.^2+u3.^2);
        uf= u1./rf;
        vf= u2./rf;
        if PFLAG %If user input custom pupil
            %             cosmaxp=real(-qb./(2*qxb.*rh).*(sqrt(rh.^2+qb.^2/4).*(ca+cb)+abs(sqrt(rh.^2+qb.^2/4).*(ca-cb)+sgn.*q3b)));
            %             phixp=real(acos(complex(cosmaxp)));
            %             phixminp= phixp;
            %             phixmaxp= mod(2*pi-phixp,2*pi);
            %             phip=[0.5:Nph-.5]./Nph.*(phixmaxp-phixminp)+phixminp;
            %             s1_p= rh.*cos(phip);
            %             s2_p= rh.*sin(phip);
            cosmaxp=real(-qb./(2*qxb.*rh).*(sqrt(rh.^2+qb.^2/4).*(ca+cb)+abs(sqrt(rh.^2+qb.^2/4).*(ca-cb)+sgn.*q3b)));
            phixp=real(acos(complex(cosmaxp)));
            phixminp= phixp;
            phixmaxp= mod(2*pi-phixp,2*pi);
            phip=[0.5:Nph-.5]./Nph.*(phixmaxp-phixminp)+phixminp;
            s1_p= rh.*cos(phip);
            s2_p= rh.*sin(phip);
            u1_p= -sgn.*q1b/2-q2b.*s2_p./qxb+q1b.*q3b.*s1_p./(qxb.*qb);
            u2_p= -sgn.*q2b/2+q1b.*s2_p./qxb+q2b.*q3b.*s1_p./(qxb.*qb);
            u3_p= -sgn.*q3b/2-qxb.*s1_p./qb;
            rf_p= sqrt(u1_p.^2+u2_p.^2+u3_p.^2);
            uf_p= u1_p./rf_p;
            vf_p= u2_p./rf_p;
            if GPU
                sample=interpn(Uf,Vf,Rf,Sf,uf,vf,rf,'linear',0)./rf.^2;
                pupil1=interpn(Uf,Vf,Rf,Pf,uf_p,vf_p,rf_p,'linear',0);
                if BORN
                    pupil2=interpn(Uf,Vf,Rf,Pf,uf,vf,rf,'linear',0);
                    if PCX
                        pupil1_ang=interpn(Uf,Vf,Rf,Pf_ang,uf_p,vf_p,rf_p,'linear',0);
                        pupil2_ang=interpn(Uf,Vf,Rf,Pf_ang,uf,vf,rf,'linear',0);
                    end
                end
            else
                sample=Fsrc(uf,vf,rf)./rf.^2;
                pupil1=Fppl(uf_p,vf_p,rf_p);
                if BORN
                    pupil2=Fppl(uf,vf,rf);
                    if PCX
                        pupil1_ang=Fppl_ang(uf_p,vf_p,rf_p);
                        pupil2_ang=Fppl_ang(uf,vf,rf);
                    end
                end
                pupil1(isnan(pupil1))=0;
                pupil2(isnan(pupil2))=0;
                sample(isnan(sample))=0;
            end
            if BORN
                if PCX %complex pupil
                    fun= sample.*(pupil1.*exp(1i.*pupil1_ang)).*(pupil2.*exp(1i*pupil2_ang))./qb;
                else
                    fun= sample.*(pupil1).*(pupil2)./qb;
                end
            else %RYTOV
                fun= sample.*pupil1.^2./qb;
            end
        else
            if GPU %BORN and RYTOV same
                sample=interpn(Uf,Vf,Rf,Sf,uf,vf,rf,'linear',0)./rf.^2;
            else
                sample=Fsrc(uf,vf,rf)./rf.^2;
                sample(isnan(sample))=0;
            end
            fun= sample./qb;
        end
        drh=(rhmax-rhmin)./Nr;
        dr=(vmax-vmin)./Nr;
        dph=(phixmax-phixmin)./Nph;
        tint=trapz(trapz(real(fun.*rh.*drh.*dph./pi),1),2); %check constant here
        tab(bInd)=tab(bInd)+shiftdim(tint,2);
        tph(bInd)=tph(bInd)+shiftdim(sgn.*tint,2);
    end
    if mod(bb,50)==0
        disp([num2str(bb),'/',num2str(Nbatch)])
        toc(batchTime)
    end
    if GPU
        Tab(mInd)=gather(tab(:));
        Tph(mInd)=gather(tph(:));
    else
        Tab(mInd)=tab(:);
        Tph(mInd)=tph(:);
    end
end
if any(strcmpi(options,'verbose'))
    toc(myFastTime);
end
TF=Tab+1i*Tph;
if any(strcmpi(options,'plot'))
    if any(strcmpi(options,'PSF'))
        figure(1); subplot(1,2,1); imagesc(z,x,squeeze(real(TF(:,floor(Ny/2),:)))); axis equal tight; colorbar;  title(['Absorption, NA=',num2str(NAp)]);
        subplot(1,2,2); imagesc(z,x,squeeze(imag(TF(:,floor(Ny/2),:)))); axis equal tight; colorbar; title(['Phase']);
        set(gcf,'Position',[1          31        1920         973]); drawnow;
    else
        figure(1); subplot(1,2,1); imagesc(squeeze(wl),ul,squeeze(real(TF(:,floor(Ny/2),:)))); axis equal tight; colorbar;  title(['Absorption, NA=',num2str(NAp)]);
        subplot(1,2,2); imagesc(squeeze(wl),ul,squeeze(imag(TF(:,floor(Ny/2),:)))); axis equal tight; colorbar; title(['Phase']);
        set(gcf,'Position',[1          31        1920         973]); drawnow;
    end
end
if any(strcmpi(options,'gridOut'))
    varargout{1}=Ui;
    varargout{2}=Vi;
    varargout{3}=Wi;
end
if any(strcmpi(options,'upsample'))
    Fta=griddedInterpolant(Ul3,Vl3,Wl3,Tab,'linear','none');
    Ftp=griddedInterpolant(Ul3,Vl3,Wl3,Tph,'linear','none');
    Ta=Fta(Ui,Vi,Wi); Ta(isnan(Ta))=0;
    Tp=Ftp(Ui,Vi,Wi); Tp(isnan(Tp))=0;
    TF=Ta+1i*Tp;
end

if any(strcmpi(options,'PSF'))
    ta=real(ifftshift(ifftn(ifftshift(Ta))));
    tp=imag(ifftshift(ifftn(ifftshift(Tp))));
    TF=ta+1i*tp;
end
end

