function I = loadOBMzStack_indicateEndian(dataPath,J, label,Nz,dims,centers,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Patrick Ledwig 12/12/19 Georgia Institute of Technology
% Zhenmin - Li 12/12/22 
% loadOBM: loads binary obm images taken with BlinkMoveXY.vi, for when you
% have to do the averaging yourself
% INPUTS:
%  path - location of data
%  J - struct that contains information of bin files in running folder
%  label - label of data, or index of label
%  Nz - number of slices to load
%  [Nx,Ny] - size of the ROI
%  [Cx,CY] - center pixels of ROI
%  options - optional specifications:
%   'visible': show dpc images as they are loaded
%  order - optional non-default loading order e.g. [1,3,4,2]. [1 2 3 4] is
%   default
% OUTPUTS:
%  I - cell containing images (I{color}{side})
%  J - structure containing information for the images
% EXAMPLE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=varargin{1};
if length(varargin)>=2
    order=varargin{2};
else
    order=[1,2,3,4];
end
if length(varargin)>=3
    which=varargin{3};
else
    which=1:Nz;
end
if length(varargin)>=4
    bigEndian = varargin{4};
end

Npix = dims(1);
Nx=dims(1); Ny=dims(2);
Cx=centers(1); Cy=centers(2);
ROIx=[Cx-ceil(Nx/2)+1:Cx+floor(Nx/2)];
ROIy=[Cy-ceil(Ny/2)+1:Cy+floor(Ny/2)];

% j=1;
% files=dir(dataPath);
% files= files(3:end);
% [~,ndx,~] =  natsort({files.name});
% files = files(ndx);


uID=unique({J.ID});
% uN=unique([J.N]);
uX=unique([J.X]);
uY=unique([J.Y]);
uZ=flip(sort(unique([J.Z])));
uZ=uZ(~isnan(uZ)&(uZ~=0));
%%
if ~isnumeric(label)
    J=J(contains({J.ID},label));
else
    J=J(strcmp({J.ID},uID{label}));
end
tic;
uX=uX(~isnan(uX));
uY=uY(~isnan(uY));

% 
% for i=1:length(uX)
%     x=uX(i);
%     for j=1:length(uY)
        % y = uY(j);

pp=1;
for p=which
    for c=1:4
        Ic{c}=0;
    end


    for c=1:4
        m=4*(p-1)+c;
        fid=fopen([dataPath,filesep,J(m).fname]);

        if ~J(m).mlb
            I0=fread(fid,[Npix,Npix],'uint16',0);
        else
            I0=fread(fid,[Npix,Npix],'uint16',0,'b');
        end
        
        I0(1:5,1)=I0(1:5,2);
        varI0=var(I0(:));
        I0=I0(ROIx,ROIy);
        Ic{c}=I0+I0./varI0;
        fclose(fid);
    end
    IL_R=Ic{order(1)};
    IR_R=Ic{order(2)};
    IL_G=Ic{order(3)};
    IR_G=Ic{order(4)};
    if any(strcmpi(options,'visible'))
        figure(11); imagesc((IL_R-IR_R)./(IL_R+IR_R)); axis equal tight; colorbar; title(num2str(p)); drawnow;
        %%%% PALOMA DEBUG %%%%%
        figure(2); subplot(221); imagesc(IL_R); subplot(222); imagesc(IL_G); subplot(223); imagesc(IR_R); subplot(224); imagesc(IR_G)
        colormap bone; title(num2str(p)); drawnow;
        pause(.01)
        %%%%%%%%%%%%%%%%%%%%%%%%%
    end

    
    I{1}{1}(:,:,pp)=IL_R;
    I{1}{2}(:,:,pp)=IR_R;
    I{2}{1}(:,:,pp)=IL_G;
    I{2}{2}(:,:,pp)=IR_G;
    pp=pp+1;
    clc;
    disp([num2str(p),'/',num2str(max(which))]);
    toc;
end
%     end
% end



