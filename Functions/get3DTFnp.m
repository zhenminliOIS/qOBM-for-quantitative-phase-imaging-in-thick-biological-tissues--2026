function [TF,B,varargout]=get3DTFnp(X,Y,Z,S,NA,lambda,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Ledwig 12/4/2019 Georgia Institute of Technology
% get3DTFup(X,Y,Z,S,P,options) outputs a transfer function for a microscope
% with a specified angular source distribution and exit pupil diameter for
% a quasi-monochromatic illumination at a specified wavelength
% DEPENDENCIES: 
%  parallel computing toolbox (optional) - for GPU capability 
%  cspad.m - pads the pupil for circle shifting
% INPUTS:
% X,Y,Z: 3D coordinate grid in ndgrid format in distance, e.g. [um]
% S: source distribution either as a 2D grid of values (0-1) in normalized angle
%  space (i.e.u=linspace(-1,1,Nu), v=linspace(-1,1,Nv)) or a single valued
%  number between 0 and 1 representing the NA of an entrance pupil with
%  even illumination
% lambda - wavelength to be used in same unit as X,Y,Z e.g. [um]
% NA: single valued number between 0 and 1 representing the exit pupil NA
%  with uniform transmission
% options: strings representing options to vary operation of the function
%   'PSF' - outputs instead the PSF of the microscope
%   'verbose' - outputs periodic timing updates on the completion
%    percentage of the calculation
%   'GPU' - perform computations with a gpu (requires parallel computing
%    toolbox). this speeds up computation considerably for larger grids
%   'plot' - optional plot showing a cross section through y
%   'gridOut' - outputs the non-normalized frequency space grid
%   'upsample' - lets you produce a 3D transfer function with a more highly
%   sampled grid, and then downsample it to produce a smoother TF
% updims: (optional) upsampling dimensions
% OUTPUTS:
%  TF - The 3D transfer function for the microscope. The absorption TF to be
%  given in the real part and the phase TF to be given in the imaginary
%  part
% EXAMPLE:
%  clear all; close all; clc;
%  Nx=127; Ny=127; Nz=63;
%  dx=630e-3/4;
%  dy=630e-3/4;
%  dz=630e-3/2;
%  x=((0:(Nx-1))-floor(Nx/2))*dx;
%  y=((0:(Ny-1))-floor(Ny/2))*dy;
%  z=((0:(Nz-1))-floor(Nz/2))*dz;
%  [X,Y,Z]=ndgrid(x,y,z);
%  uMax=1/dx;
%  vMax=1/dy;
%  u=linspace(-1,1,Nx); v=linspace(-1,1,Ny);
%  lambda=.630;
%  [U,V]=ndgrid(u,v);
%  S=double(sqrt(U.^2+V.^2)<=.65);
%  S=S.*(U>0)+S.*(U<=0);
%  % S=.4;
%  NA=.65;
%  lambda=.630;
%  options={'plot','verbose'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
if any(strcmpi(options,'verbose'))
    VERBOSE=1;
else
    VERBOSE=0;
end


k=2*pi/lambda;
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
if any(strcmpi(options,'upsample'))
    Nx0=updims(1); Ny0=updims(2); Nz0=updims(3);
    x0=((0:(Nx0-1))-floor(Nx0/2))*dx;
    y0=((0:(Ny0-1))-floor(Ny0/2))*dy;
    z0=((0:(Nz0-1))-floor(Nz0/2))*dz;
    
    %Define un-normalized 2d angle space for upsampling
    u=((0:(Nx0-1))-floor(Nx0/2))*uMax/Nx0;
    v=((0:(Ny0-1))-floor(Ny0/2))*vMax/Ny0;
    w=((0:(Nz0-1))-floor(Nz0/2))*wMax/Nz0;
    
    %Define the down-sample space for the output
    ui=((0:(Nx-1))-floor(Nx/2))*uMax/Nx*lambda;
    vi=((0:(Ny-1))-floor(Ny/2))*vMax/Ny*lambda;
    wi=((0:(Nz-1))-floor(Nz/2))*wMax/Nz*lambda;
    [Ui,Vi,Wi]=ndgrid(ui,vi,wi);


    Nx=Nx0; Ny=Ny0; Nz=Nz0;
else
    u=((0:(Nx-1))-floor(Nx/2))*uMax/Nx;
    v=((0:(Ny-1))-floor(Ny/2))*vMax/Ny;
    w=((0:(Nz-1))-floor(Nz/2))*wMax/Nz;
    
    %Define the down-sample space for the output
    ui=((0:(Nx-1))-floor(Nx/2))*uMax/Nx*lambda;
    vi=((0:(Ny-1))-floor(Ny/2))*vMax/Ny*lambda;
    wi=((0:(Nz-1))-floor(Nz/2))*wMax/Nz*lambda;
    [Ui,Vi,Wi]=ndgrid(ui,vi,wi);
end



%Define normalized 2D angular space from image space
ul=u*lambda; vl=v*lambda;
dul=ul(2)-ul(1); dvl=vl(2)-vl(1);
[Ul,Vl]=ndgrid(ul,vl);
%Define 3rd dimension for 3D angle space
wl=w*lambda; dwl=wl(2)-wl(1);
[Ul3,Vl3,Wl3]=ndgrid(ul,vl,wl); 

%For an arbitrary source distribution function
if numel(S)~=1
    %Resample entrance pupil into normalized spatial frequency units
    u=linspace(-1,1,size(S,1)); v=linspace(-1,1,size(S,2));
    [U,V]=ndgrid(u,v);
    Fs=griddedInterpolant(U,V,S,'linear','none'); % frequency Source 
    S=Fs(Ul,Vl); S(isnan(S))=0; 
%     S=S./(sum(S(:)));  
    
    %Sort out entrance and exit pupil
    R=sqrt(Ul.^2+Vl.^2);
    Pi=single(R<=NA);
    [~,Rind]=sort(R(:));
    P0=single(R<=(NA+1));
    B=sum(S(:).*Pi(:)).*dul.*dvl;
    
    %Parameters to make to computation more concise
    indP=(find(P0));
    indPSort=Rind(1:numel(indP));
    MidU=floor(Nx/2);
    MidV=floor(Ny/2);
    
    %Initialize varaibles to be used
    Ta=zeros(Nx,Ny,Nz); %absorption transfer function
    Tp=zeros(Nx,Ny,Nz); %phase transfer function
    wl=reshape(wl,[1,1,length(wl)]);
    if GPUFLAG
        Ta=gpuArray(Ta);
        Tp=gpuArray(Tp);
        Pi=gpuArray(Pi);
        S=gpuArray(S);
        Ul=gpuArray(Ul);
        Vl=gpuArray(Vl);
        wl=gpuArray(wl);
        Ul3=gpuArray(Ul3);
        Vl3=gpuArray(Vl3);
    end
%     if VERBOSE
%         tfTime=tic;
%     end
    %this code loops over each angular frequency in the exit pupil space
    %and computes the depth transfer function as given by Streibl 1985 for
    %an arbitrary source distribution. The delta function in the above is
    %represented as a unit step in normalized w (wave vector component 3)
%     f_waitbar = waitbar(0,'get3DTFnp, Please wait...');
%     pause(.5)  
    for nn=1:numel(indP)
        mm=indPSort(nn);
        [um,vm]=ind2sub([Nx,Ny],mm);
        uu=um-MidU-1;vv=vm-MidV-1;
        u1=-ceil(uu/2);v1=-ceil(vv/2);
        u2=floor(uu/2);v2=floor(vv/2);
        P1=circshift(Pi,[u1,v1]);
        P1=cspad(P1,[u1,v1]);
        P2=circshift(Pi,[u2,v2]);
        P2=cspad(P2,[u2,v2]);
        S1=circshift(S,[u1,v1]);
        S1=cspad(S1,[u1,v1]);
        S2=circshift(conj(S),[u2,v2]);
        S2=cspad(S2,[u2,v2]);
        FA=1./(4*pi*dwl)*P1.*P2.*(S1+S2);
        FP=1./(4*pi*dwl)*P1.*P2.*(S1-S2);
        wsum1=1-(Ul3+.5*Ul(um,vm)).^2-(Vl3+.5*Vl(um,vm)).^2;
        wsum2=1-(Ul3-.5*Ul(um,vm)).^2-(Vl3-.5*Vl(um,vm)).^2;
        wu1=sqrt(wsum1.*(wsum1>=0));
        wu2=sqrt(wsum2.*(wsum2>=0));
        w1=-(wu1-wu2)>=-(wl+dwl/2);
        w2=-(wu1-wu2)<-(wl-dwl/2);
        wM=w1.*w2;
        Ta(um,vm,:)=sum(sum(wM.*FA,1),2)*dul*dvl;
        Tp(um,vm,:)=sum(sum(wM.*FP,1),2)*dul*dvl;

        if (mod(nn,1000)==0)&&(VERBOSE)
%             tftime = toc(tfTime);
%             disp([num2str(nn),'/',num2str(numel(indP))]);
%             fprintf(repmat('\b',1,count));
            
%             waitbar(nn/numel(indP),f_waitbar, [num2str(nn),'/',num2str(numel(indP))]);
%             pause(0.5)
%             tfTime=tic;
        end
    end
%     delete(f_waitbar);

    if GPUFLAG
        Ta=gather(Ta);
        Tp=gather(Tp);
        Pi=gather(Pi);
        S=gather(S);
        Ul=gather(Ul);
        Vl=gather(Vl);
        wl=gather(wl);
        Ul3=gather(Ul3);
        Vl3=gather(Vl3);
    end
    TF=Ta+1i*Tp;
    %For a circular entrance pupil  w/ uniform illumination, from Streibl 1985
elseif numel(S)==1
    NAs=S;
    ps=NAs;
    pp=NA;
    p=sqrt(Ul3.^2+Vl3.^2);
    n=Wl3.*lambda;
    Ta=1/(2*pi*p).*real((.5*(pp^2+ps^2)-...
        1/4*p.^2-(n./(lambda*p)).^2-...
        abs(n/lambda-.5*(pp^2-ps^2))).^.5+...
        (.5*(pp^2+ps^2)-...
        1/4*p.^2-(n./(lambda*p)).^2-...
        abs(n/lambda+.5*(pp^2-ps^2))).^.5);
    Tp=1./(2*pi*p).*real((.5*(pp^2+ps^2)-...
        1/4*p.^2-(n./(lambda*p)).^2-...
        abs(n/lambda-.5*(pp^2-ps^2))).^.5-...
        (.5*(pp^2+ps^2)-...
        1/4*p.^2-(n./(lambda*p)).^2-...
        abs(n/lambda+.5*(pp^2-ps^2))).^.5);
    
    Ta(isnan(Ta))=0;
    Tp(isnan(Tp))=0;
    TF=Ta+1i*Tp;
    B=pi.*min(NA,NAs).^2;
end
if any(strcmpi(options,'plot'))
    if any(strcmpi(options,'PSF'))
        figure(1); subplot(1,2,1); imagesc(z,x,squeeze(real(TF(:,floor(Ny/2),:))));  colorbar;  title(['Absorption, NA=',num2str(NA)]);
        subplot(1,2,2); imagesc(z,x,squeeze(imag(TF(:,floor(Ny/2),:))));  colorbar; title(['Phase']);
        set(gcf,'Position',[1          31        1920         973]); drawnow;
    else
        figure(1); subplot(1,2,1); imagesc(squeeze(wl),ul,squeeze(real(TF(:,floor(Ny/2),:))));  colorbar;  title(['Absorption, NA=',num2str(NA)]);
        subplot(1,2,2); imagesc(squeeze(wl),ul,squeeze(imag(TF(:,floor(Ny/2),:))));  colorbar; title(['Phase']);
        set(gcf,'Position',[1          31        1920         973]); drawnow;
    end
end
if any(strcmpi(options,'gridOut'))
    varargout{1}=Ui;
    varargout{2}=Vi;
    varargout{3}=Wi;
end
if any(strcmpi(options,'upsample'))
    Fta=griddedInterpolant(Ul3,Vl3,Wl3,Ta,'linear','none');
    Ftp=griddedInterpolant(Ul3,Vl3,Wl3,Tp,'linear','none');
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