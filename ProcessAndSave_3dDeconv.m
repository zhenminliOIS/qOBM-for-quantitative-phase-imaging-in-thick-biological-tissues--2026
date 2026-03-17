clc
close all
clear

%% Set Paths
Path = 'Data\Raw Human Brain';
OutPath = 'Y:\robles\OBM\qOBM protocol\Data\3D deconv divide 2B\';
mcDir = 'MC\3D\';

%% Processing Settings
simulateTF = 1;
gamma = -1.8000;
Kx = 2;
Ky = 2;
Ox = 256;
Oy = 256;
GPU = 1;
WD = 2;
filter = 0;
n0 = 1.4; % average of RI

WVL = [720 720];
mag = 60;
type = 'Brain';
inputMag = 60;
invertL = 0;
invertR = 0;
defTH = 45;
runOnlyCheck = 0;
runOnly = '2021_Case10_SideA_60x_45deg_Z4';
NoAngles = 0;
bigEndian = 1;

%% Add Function Paths and define dimensions
addpath('Functions')
dims = [2048,2048];
paddingLayer = 1;

l1 = 10; 
h1 = 11;     % parameters for paper diffuse           
crop = 15;

lambda = WVL(1) / 1000;   % wavelength in um

%% Set up saving paths
dataSavePath = fullfile(OutPath, 'saved_data');
tfSavePath   = fullfile(OutPath, 'saved_tf');
if ~exist(dataSavePath, 'dir')
    mkdir(dataSavePath)
end
if ~exist(tfSavePath, 'dir')
    mkdir(tfSavePath)
end

%% >>>>>>>>>> Start of 3D deconvolution <<<<<<<<<<<<
% Check for subfolders in the given directory
thisFolder = false;
files = dir(Path);
dirFlags = [files.isdir] & ~strcmp('.', {files.name}) & ~strcmp('..', {files.name});
Folders = files(dirFlags);

if isempty(Folders)
    Folders = 0;
    thisFolder = true;
end

%% Loop over each folder
for d = 1:length(Folders)
    if thisFolder
        dataPath = Path;
        thisPath = Path;
        files = dir(dataPath);
    else
        dataPath = fullfile(Path, Folders(d).name);
        thisPath = dataPath;
        files = dir(dataPath);
    end

    % Read files and determine unique labels for processing
    [Jfiles, ~] = readFiles(files, runOnly, runOnlyCheck, bigEndian);
    dataLabels = unique({Jfiles.ID});

    %% Process each individual z-scan
    for pp = 1:length(dataLabels)
        thisLabel = dataLabels{pp};
        disp(['************  Processing data number ', num2str(pp), ', label = ', thisLabel, '  ***************'])
        
        filePattern = fullfile(dataPath, ['*', thisLabel, '*.bin']);
        filesLabel = dir(filePattern);
        numFiles = numel(filesLabel);
        disp(['numFiles = ', num2str(numFiles)])
        
        %% Check if this is a z-scan; if not, skip this label
        XYZmat = zeros(numFiles, 3); % matrix for x,y,z info
        indices = contains({Jfiles.ID}, thisLabel);
        thisJ = Jfiles(indices);
        for ppp = 1:numFiles
            XYZmat(ppp,1) = thisJ(ppp).X; 
            XYZmat(ppp,2) = thisJ(ppp).Y; 
            XYZmat(ppp,3) = thisJ(ppp).Z; 
        end
        % Use every fourth file as representative
        XYZmat = XYZmat(1:4:numFiles, :);
        XYZmean = mean(XYZmat, 1);
        
        if XYZmean(1) ~= XYZmat(1,1) || XYZmean(2) ~= XYZmat(1,2)
            disp('Not a z scan')
            continue  % skip non z-stack data
        end

        dz_stage = abs(round(mean(unique(diff(XYZmat(:,3)))), 3, 'significant') * 1000);
        if dz_stage == 0
            dz_stage = 2;
        end

        %% >>>>>>>>>> Chunk Processing <<<<<<<<<<<<
        [mag, NA, cmin, cmax, Type, ~, ~, theta, LD, freq] = getFileInformation(thisLabel, type, mag, defTH);
        TH = theta;
        dx = .138 * (60 / mag);                                 
        dy = .138 * (60 / mag); 
        
        % Correction factor based on NA
        switch NA
            case 0.7
                correction = 1.58;
            case 0.6
                correction = 1.45;
            case 0.45
                correction = 1.515;
            case 0.8
                correction = 1.7;
            otherwise
                correction = 1;
        end

        dz = abs(correction * dz_stage); 
        dvol = dx * dy * dz;
        
        % Determine number of images and imaging dimensions
        N_img = floor(numel(filesLabel) / 4);  
        Mx = floor((dims(1) - (Kx + 1) * Ox) / Kx);
        My = floor((dims(2) - (Ky + 1) * Oy) / Ky);
        Nx = Mx + 2 * Ox;
        Ny = My + 2 * Oy;
        
        % Define real-space axes
        x = ((0:(Nx-1)) - floor(Nx/2)) * dx;
        y = ((0:(Ny-1)) - floor(Ny/2)) * dy;
        z = ((0:(N_img-1)) - floor(N_img/2)) * dz;
        [X, Y, Z] = ndgrid(x, y, z);
        
        % Frequency coordinates for source distribution
        u0 = linspace(-1, 1, Nx);
        v0 = linspace(-1, 1, Ny);
        [U0, V0] = ndgrid(u0, v0);
        p0 = sqrt(U0.^2 + V0.^2);
        
        %% Get source distribution S based on sample type
        switch type
            case 'Brain'
                loadLbl = ['BRN_WD', num2str(WD), '_', num2str(WVL(1)), 'nm_th', num2str(TH), '_sep', num2str(LD), 'mm'];
                tfFile = fullfile(mcDir, [loadLbl, '_Kx', num2str(Kx), '_Ky', num2str(Ky), '_trunked.mat']);
                if exist(tfFile, 'file')
                    load(tfFile, 'S');
                else
                    S = brainDiffuse([Nx, Ny], mcDir, loadLbl);
                    save(tfFile, 'S');
                end
            case 'Paper'
                S = paperDiffuse([Nx, Ny], l1, h1, theta * pi/180);
        end

        % Resize and normalize S
        S = imresize(S, length(x) ./ size(S, 1));
        S0{1}{1} = S;
        S0{1}{1}(isnan(S0{1}{1})) = 0;
        S0{1}{1} = S0{1}{1} ./ sum(S0{1}{1}(:));
        
        % Create flipped and rotated versions of S
        S0{1}{2} = flip(S0{1}{1}, 1);
        S0{2}{1} = rot90(S0{1}{1});
        S0{2}{2} = flip(rot90(S0{1}{1}), 2);
        
        %% >>>>>>> 3D Microscope TF Calculation <<<<<<<<
        if ~exist('TFdpc3D', 'var')
            if GPU
                opts = {'GPU','verbose','gridOut','upsample'};
            else
                opts = {'verbose','plot','gridOut'};
            end

            [TF3d1, B, Ul, Vl, Wl] = get3DTFnpbb(X, Y, Z, S0{1}{1}, NA, lambda, opts, [255 255 size(X,3)]);  
            ul = sort(unique(Ul)); 
            dul = ul(2) - ul(1);
            vl = sort(unique(Vl)); 
            dvl = vl(2) - vl(1);
            wl = sort(unique(Wl)); 
            dwl = wl(2) - wl(1);
            
            if GPU
                opts = {'GPU','verbose','gridOut','upsample'};
            else
                opts = {'verbose'};
            end
            [TF3d2, ~] = get3DTFnpbb(X, Y, Z, S0{1}{2}, NA, lambda, opts, [255 255 size(X,3)]); 
            TFdpc3D = (imag(TF3d1) - imag(TF3d2))./ (2 * B);
            
            
            save(fullfile(tfSavePath, [thisLabel, '_TFdpc3D.mat']), 'TFdpc3D', 'S')
            clear TF3d2 TF3d1;
        end
        disp('Load 3D TF finished')
        
        %% Load stored data
        disp('Loading the data')
        I3full = loadOBMzStack_indicateEndian(dataPath, Jfiles, thisLabel, N_img, dims, [1024,1024], {'non visible'}, [1,2,3,4], 1:N_img);
        
        %% Process data chunk by chunk
        for c_x = 1:Kx
            for c_y = 1:Ky
                disp(['cx = ', num2str(c_x), '  cy = ', num2str(c_y)]);
                chunkLabel = [thisLabel, '_cx', num2str(c_x), '_cy', num2str(c_y)];
                
                % Determine chunk centers
                Cx = floor(Mx/2) + c_x * (Ox) + (c_x - 1) * Mx; 
                Cy = floor(My/2) + c_y * (Oy) + (c_y - 1) * My;
                
                dz3 = 1;
                wMax3 = 1 / dz3;
                ul3 = sort(unique(Ul)); 
                vl3 = sort(unique(Vl));
                wl3 = ((0:(N_img-1)) - floor(N_img/2)) * wMax3 / N_img * lambda;
                [Ul3, Vl3, Wl3] = ndgrid(ul3, vl3, wl3);
                
                indexMat = reshape(1:prod([2048,2048,N_img]), [2048,2048,N_img]);
                chunkInds = indexMat(Cx - floor(Nx/2) + 1 : Cx + floor(Nx/2), ...
                                     Cy - floor(Ny/2) + 1 : Cy + floor(Ny/2), 1:N_img);
                clear indexMat;
                
                % Chunk the images for each channel and side
                for col = 1:2
                    for ii = 1:2
                        I3c{col}{ii} = I3full{col}{ii}(chunkInds);
                    end
                end
                clear chunkInds;
                
                % Calculate DPC images for the chunk
                Idpcs{1} = (I3c{1}{1} - I3c{1}{2}) ./ (I3c{1}{1} + I3c{1}{2});
                Idpcs{2} = (I3c{2}{1} - I3c{2}{2}) ./ (I3c{2}{1} + I3c{2}{2});
                
                if filter
                    [filtH] = RadialFiltGen(Nx, Ny, -2, 15, 6);
                    [filtL] = RadialFiltGen(Nx, Ny, 2, Nx/3, 9);
                    Filt = repmat(filtH .* filtL, 1, 1, N_img);
                    for col = 1:2
                        Idpcs{col} = real(ifft2(ifftshift(Filt, 1, 2) .* fft2(Idpcs{col})));
                    end
                end
                
                %% 1st: Get or calculate incident angle (TH)
                disp('Get or calculate TH')
                if ~NoAngles  % use stored TH
                    fid = fopen(fullfile(thisPath, [thisLabel, '_angles.txt']));
                    TH = num2cell(fscanf(fid, '%f'));
                    fclose(fid);
                    % Adjust the TH for each channel
                    TH{1} = TH{1} + pi;
                    TH{2} = TH{2} + pi;
                    
                    % Plot to verify incident light direction
                    u0 = 0; v0 = 0;
                    u1 = 0; 
                    v1 = floor(Nx/2) * dul;
                    figure(29);
                    for col = 1:2   
                        ur = u1 * cos(TH{col} - pi/2) + v1 * sin(TH{col} - pi/2); 
                        vr = u1 * -sin(TH{col} - pi/2) + v1 * cos(TH{col} - pi/2);
                        subplot(1,2,col);
                        imagesc(ul, vl, squeeze(log(abs(sum(fftshift(fft2(Idpcs{col}),1), 3)))));
                        axis equal tight xy;
                        colormap gray;
                        hold on;
                        plot([u0, ur], [v0, vr], '-r');
                        hold off;
                    end
                else
                    for col = 1:2
                        for n = 1:N_img
                            FF(:,:,n) = fftshift(fft2(Idpcs{col}(:,:,n)));
                        end
                        ffsum{col} = sum(FF,3);
                        clear FF;
                        IDff = ifft2(ifftshift(ffsum{col}));
                        TH{col} = ShearAngle(IDff);
                    end
                end
                
                %% 2nd: Rotate the transfer function using TH
                for col = 1:2
                    Ftdpc = griddedInterpolant(Ul, Vl, Wl, TFdpc3D);
                    Ulr = cos(TH{col}) .* Ul3 + sin(TH{col}) .* Vl3;
                    Vlr = -sin(TH{col}) .* Ul3 + cos(TH{col}) .* Vl3;
                    TD3d_Rot{col} = Ftdpc(Ulr, Vlr, Wl3);
                    clear Ftdpc
                end
                
                %% DPC Reconstruction
                if invertL
                    TD3d_Rot{1} = -TD3d_Rot{1};
                end
                if invertR
                    TD3d_Rot{2} = -TD3d_Rot{2};
                end

                B = max(abs(TD3d_Rot{1}(:))) / 0.043;
                TD3d_Rot{1} = TD3d_Rot{1} ./ B .*2;
                TD3d_Rot{2} = TD3d_Rot{2} ./ B .*2;
                
                save(fullfile(tfSavePath, [thisLabel, '_TD3d_Rot.mat']), 'TD3d_Rot')
                


                alpha = 10^(gamma);
                for col = 1:2                       
                    FD{col} = fftshift(fftn(Idpcs{col}));
                end

                Jnum = zeros(size(Ul3));
                Jden = zeros(size(Ul3));
                
                for col = 1:2
                    TD3d_Rot{col}(isnan(TD3d_Rot{col})) = 0;
                
                    Jnum = Jnum + ...
                        conj(-1i * FD{col}) .* 2 .* TD3d_Rot{col} .* ((-1)^(col+1));
                
                    Jden = Jden + ...
                        abs(2 * TD3d_Rot{col}).^2 + alpha;
                end
                
                J = conj(Jnum) ./ Jden;
                j = ifftn(ifftshift(J)) / dvol;
                
                k = 2 * pi * n0 / lambda;
                n = real(sqrt(-j/k^2 + n0^2));
                phaseVol3D = 2 * pi / lambda * (n - n0)*dz;
                
                save(fullfile(dataSavePath, [chunkLabel, '.mat']), 'phaseVol3D')
            end % for c_y
        end % for c_x

        %% Clear temporary variables for current label
        clear I3full I3c Tpr Idpcs TD3d_Rot TFdpc3D

        %% Image and Video Export
        outPathImg = fullfile(OutPath, 'Images');
        outPathMov = OutPath;
        if ~exist(outPathImg, 'dir')
            mkdir(outPathImg);
        end

        % Recall parameters for image reconstruction
        Mx = floor((2048 - (Kx + 1) * Ox) / Kx);  
        My = floor((2048 - (Ky + 1) * Oy) / Ky);
        Nx = Mx + Ox;   
        Ny = My + Oy;
        NNx = 2048 - Ox; 
        NNy = 2048 - Oy; 

        % Load and assemble chunks into phaseIMs
        phaseIMs = zeros(NNx, NNy, N_img);
        for c_x = 1:Kx
            for c_y = 1:Ky
                Cx = floor(Mx/2) + floor(Ox/2) + (c_x - 1) * (Mx + Ox);
                Cy = floor(My/2) + floor(Oy/2) + (c_y - 1) * (My + Oy);
                chunkLabel = [thisLabel, '_cx', num2str(c_x), '_cy', num2str(c_y)];
                load(fullfile(dataSavePath, [chunkLabel, '.mat'])); 
                phaseIMs(Cx - floor(Nx/2) + 1 : Cx + floor(Nx/2), ...
                         Cy - floor(Ny/2) + 1 : Cy + floor(Ny/2), :) = ...
                         phaseVol3D(floor(Ox/2) + 1 : Mx + floor(1.5 * Ox), ...
                                    floor(Oy/2) + 1 : My + floor(1.5 * Oy), :);
                delete(fullfile(dataSavePath, [chunkLabel, '.mat']));
            end
        end
        
        phaseIMs = phaseIMs(crop:end-crop+1, crop:end-crop+1, :);
        
        %% Save Image Frames
        [cmin, cmax] = percentileScaleQuick(phaseIMs(:), .01, .999);
        for nn = 1:N_img 
            imrgb = phaseIMs(:,:,nn);
            imrgb(imrgb < cmin) = cmin; 
            imrgb(imrgb > cmax) = cmax;
            mapped_image = (double(imrgb) - cmin) ./ (cmax - cmin);
            mapped_image = mapped_image * 256;
            mapped_image = uint8(mapped_image);
            mapped_image = ind2rgb(mapped_image, bone(256));
            imwrite(mapped_image, fullfile(outPathImg, [thisLabel, '_', num2str(nn), '.tif']));
        end
        
        save(fullfile(OutPath, [thisLabel, '_Full.mat']), 'phaseIMs', '-v7.3');
        clear Jnum TFdpc3D TF3dReal Ul Ul3 Uli Ulr Vl Vl3 Vli Vlr Wl Wl3 X Y Z Ftp;
    end % for pp (dataLabels)
end % for d (Folders)