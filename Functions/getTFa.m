% Paloma Casteleiro - 09.15.18
% Modification Zhenmin - 05.12.24

function H_phase = getTFa(TH,S,fitobject,NA,M,dim,WVL)

h = waitbar(0/11,'Initializing Monte Carlo Data');
for col = 1:2
    % Step 1: Initialize and Load MC Data
    waitbar((5*(col-1)+1)/11,h,'Initializing Monte Carlo Data');
    
    % Create Grid
    ui = linspace(-1,1,dim(1));
    vi = linspace(-1,1,dim(2));
    du = ui(2)-ui(1);
    dv = vi(2)-vi(1);
    [Ui,Vi] = ndgrid(ui,vi);
    
    % Step 2: Produce Transfer Function by Multiplying in the Fourier
    %         Domain Instead of Convolving to Increase Computational Speed
    waitbar((5*(col-1)+2)/11,h,'Producing Microscope Transfer Function');
    pmax = NA;
    Ri = sqrt(Ui.^2+Vi.^2);

    % Define pupil function 
    % (source distribution S was a function input)
    Pi = double(Ri<=pmax);

    % Double the size of source distribution 
    SP = zeros(size(S).*2-1);
    % Double the size of the pupil function 
    PP = zeros(size(Pi).*2-1);

    % Define S*Pi and Pi at the centers of SP and PP, respectively
    SP(ceil(size(SP,1)/4):size(S,1)+ceil(size(SP,1)/4)-1,ceil(size(SP,2)/4):size(S,2)+ceil(size(SP,2)/4)-1) = S.*Pi;
    PP(ceil(size(PP,1)/4):size(Pi,1)+ceil(size(PP,1)/4)-1,ceil(size(PP,2)/4):size(Pi,2)+ceil(size(PP,2)/4)-1) = Pi;
    
    % Get the point spread function in the phase domain for both image
    % orientations
    H1 = ifftshift(ifft2(fft2(SP).*fft2(PP))).*du.*dv; 
    H2 = flip(H1,1);
    
    % Get the DPC transfer function 
    Hph = H1-H2;

    % Normalize the image by its background
    B = sum(S(:).*Pi(:))*du*dv;
    H_phase = Hph./B;

    % Step 3: Get Shear Angle from Image
    waitbar((5*(col-1)+3)/11,h,'Calculating Shear Angle from Image');
    th = TH{col};

    % Step 4: Create Grid in Frequency space
    waitbar((5*(col-1)+4)/11,h,'Creating Frequency Space Grid');
    lambda = WVL(col)*1e-3; % Wavelength in micrometers
    dx = .138*60/M;
    dy = .138*60/M;
    Nx = dim(1);
    Ny = dim(2);
    F = 1/dx;
    G = 1/dy;
    fi = ((1:Nx)-floor(Nx/2))*F/Nx;
    gi = ((1:Ny)-floor(Ny/2))*G/Ny;
    fu = fi.*lambda;
    gu = gi.*lambda;
    [Fu,Gu] = ndgrid(fu,gu);

    % Step 5: Rotate Transfer Function to Match Shear Angle
    waitbar((5*(col-1)+5)/11,h,'Rotating Transfer Function to Match Shear Angle');
    y2 = (-floor(size(Hph,2)/2):floor(size(Hph,2)/2))*dv;
    x2 = (-floor(size(Hph,1)/2):floor(size(Hph,1)/2))*du;
    [X2,Y2] = ndgrid(x2,y2);
    % Define X and Y by the shear angle
    Xt = Fu.*cos(th)+Gu.*sin(th);
    Yt = -Fu.*sin(th)+Gu.*cos(th);
    % Rotate TF
    F = griddedInterpolant(X2,Y2,H_phase,'linear','linear');
    Hf = F(Xt,Yt);
    % Zero-out the extra parts of SP and PP
    Hf(isnan(Hf)) = 0;
    HF{col} = Hf;
end

waitbar(1,h,'Complete, Outputting Data');
H_phase = HF;
close(h)

end
