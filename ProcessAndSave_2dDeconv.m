% function ProcessAndSave(Path,OutPath,vidOptions,options,figDisplay)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paloma Casteleiro Costa - January 2021
% Norah and Zhenmin - 2025
%  Raw data is expected to be saved sequentially in the imaging order and
%  contain certain tags in the file names. These tags may be editted in
%  line
%   - magnification (example: 60x)
%   - Image number
%   - Angle of illumination
%   - Time tag - to ensure the images are saved in the right order
%   - Z position
%   - Image name example: Sample_60x_45deg_1__Z0.004_20-11-23_14-33-34.bin
%
% INPUTS:
%   Path: input path containing raw data (.bin files)
%   OutPath: output path
%
% Options:
%   WVL - LED wavelength for each image pair
%   type - Scatterer used, intralipid is the only option in this version
%   display - Boolean
%   inputMag - Default magnification
%   invertL - Invert Transfer function of left DPC
%   invertR  - Invert Transfer function of right DPC
%   defTH - Inclination angle of fibers onto the sample
%   runOnlyCheck
%   runOnly - Run only files that contain certain characters or words
%   dontOverwrite - dont overwrite data
%   NoAngles - Boolean describing whether the shear angle has been
%                      already saved (false) or it needs to be calculated (true)
%   figDisplay: axis object to display processed images
%
% REFERENCES:
% Patrick Ledwig et. Al, Biomed. Opt. Express 10, 3605-3621 (2019)
% Patrick Ledwig et. Al, Optica 8, 6-14 (2021)
% Paloma Casteleiro Costa et. Al, Biomed. Opt. Express 12, 1621-1634 (2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Path = 'Data\Raw Human Brain\';
OutPath = 'Data\2D Deconv\';

WVL = [720 720];
mag = 60;
type = 'Brain';
display = 0;
inputMag = 60;
invertL = 1;            % invert 1st DPC TF or not 
invertR = 1;            % invert 2nd DPC TF or not
defTH = 45;             % default fiber angle
runOnlyCheck = 0;       
runOnly = '';
dontRunCheck = 0;
dontRun = '';
saveZs = 0;
dontOverwrite = 1;
NoAngles = 0;           % no saved angle, then compute illumination angle
saveMatsInd = 0;
ImOptions = 1;
InconsistentAngles = 0;
changeStartFrame = 0;
bigEndian = 1;
grab = 1;

clc

oldTH{1} = 0;
oldTH{2} = 0;
oldMag = 0;
ReloadTF = true;

addpath('Functions')
dims = [2048,2048];

% Preset paths and settings
if ~exist([OutPath,'\images\'], 'dir')
    mkdir([OutPath,'\images\'])
end
if saveMatsInd
    if ~exist([OutPath,'\IndividualMats\'], 'dir')
        mkdir([OutPath,'\IndividualMats\'])
    end
end

J = [];
saveColorbar = true;

% Detect whether there are subfolders in the input path
thisFolder = false;

files = dir(Path);
dirFlags = [files.isdir].*~strcmp('.',{files(:).name}).*~strcmp('..',{files(:).name});
Folders = files(logical(dirFlags));

if dontOverwrite
    alreadySaved = dir([OutPath,'\*.mat']);
    alreadySavedAvi = dir([OutPath,'\*.avi']);
end

if isempty(Folders)
    Folders = 0;
    thisFolder = true;
end

% Find input data that matches expected name and type
J = [];
for d = 1:length(Folders)
    if thisFolder
        thisPath = Path;
        files = dir(thisPath);
    else
        thisPath = [Path,'\',Folders(d).name];
        files = dir(thisPath);
    end

    [fIDs,J] = getfIDs(files,bigEndian);
    
    if ~NoAngles 
        fIDs = unique(fIDs);
    else
        fIDs = unique({J(:).ID});
    end

    for ii = 1:length(fIDs)
        % Prealocate bins
        phaseIM = [];
        rawIMs = [];
        allComb = [];
        positions = [];
        fID = fIDs{ii};
        % Files to fun
        Jrun = J(strcmp(fID,{J(:).ID}));
       
        % If box "Run only files containing" is checked
        if runOnlyCheck
            if ~contains(fID,runOnly,'IgnoreCase',true)
                continue;
            end
        end

        % If box "Do NOT run files containing" is checked
        if dontRunCheck
            if contains(fID,dontRun,'IgnoreCase',true)
                continue;
            end
        end

        % Set magnification and NA from the objectives available in the
        % system
        % Individual users can change getFileInformation function
        [mag, NA, cmin, cmax, Type, fitobject, S, theta, LD, freq] = getFileInformation(fID, type, mag, defTH);

        if sum(S(:))>1
            S = S./sum(S(:));
        end
       

        % Update values in Z scaled based on NA and brain RI 
        n2 = 1.4; % average RI of brain;
        zScale = (tan(asin(NA))./(tan((asin(NA))./(n2))));
        temp = num2cell([Jrun(:).Z].*zScale);
        [Jrun.Z] = temp{:};
        theseZ = unique([Jrun(:).Z]);
        ZeroZ = max(theseZ);

        if isempty(Jrun); continue; end

        %Fix issue with wrong X Y Z, (the first Z of every 4 is wrong)
        avrgs = 1;
        Jrun = table2struct(sortrows(struct2table(Jrun),{'time','N','Y','X'})); %Order struct by N
        
        % Get shear angle from pre-stored file
        if ~NoAngles 
            formatSpec = '%f';
            % Data format to read = floating-point numbers
            fid=fopen([thisPath,'\',fID,'_angles.txt']);
            % Creates a cell array of angles from _angles.txt file
            TH = num2cell(fscanf(fid,formatSpec));
            fclose(fid);

            % Check image size was saved correctly
            m = ceil((length(Jrun)/4)/2);
            if (4*(m-1)+4) > length(Jrun); continue; end
            if grab
                % This is for if the image size was correct
                fid=fopen([thisPath,'\',Jrun(4*(m-1)+1).fname]);
            else
                % This fixes an incorrectly saved image size
                fid=fopen([thisPath,'\',Jrun([Jrun.N]==(4*(m-1)+mm-1)).fname]);
            end

            if (fid==-1)
                break;
            end
            
            if Jrun(4*(m-1)+1).mlb
                % This is for correctly-sized images
                I0=fread(fid,[2048,2048],'uint16');
            else
                % This is for resized images
                I0=fread(fid,[2048,2048],'uint16',0,'b');
            end
            fclose(fid);

            % Setting minimum size of the image
            [minSize,ind] = min([dims,size(I0)]);
            % Setting maximum size of the image
            maxSize = max(size(I0));
            if ind == 3 || ind == 4; minSize = minSize - 1; end
            % Crop to smallest image size
            % dims = [minSize,minSize];
            % Or pad up to largest image size
            dims = [maxSize,maxSize];
            % Save angles on text file in raw data folder for future use,
            % or to edit

            % Skip reloading TF when all angles are the same to save on time
            if mag == oldMag && TH{1} == oldTH{1} && TH{2} == oldTH{2} && ~InconsistentAngles
                ReloadTF = false; 
            end

            if ImOptions == 1 && ReloadTF
                % Run Monte Carlo simulation for the tissue type needed. GPU expected
                if strcmp('Paper',type) 
                    % If the scatterer is paper
                    loadPaper=1;
                    H_phase = loadPaperTF(TH,mag,NA,dims,WVL,theta,LD,loadPaper);
                else
                    % If the scatterer is brain
                    H_phase = getTFa(TH,S,fitobject,NA,mag,dims,WVL);
                end
    
                % Invert boxes
                if invertL
                    H_phase{1} = -H_phase{1};
                end
                if invertR
                    H_phase{2} = -H_phase{2};
                end
            end
            ReloadTF = true;
            oldMag = mag;
            oldTH = TH;
    
            % If Save each Z individually box is checked
            % Order data by time to ensure the correct sequential order
            if ~saveZs
            % Process Z stacks first and save each Z individually
                Jrun = table2struct(sortrows(struct2table(Jrun),{'time','N','X','Y','Z','N'}));
            else
                Jrun = table2struct(sortrows(struct2table(Jrun),{'time','X','Y','Z','N'}));
            end
     
    
            % Need this dpcOnly variable for image processing and displaying
            if ~strcmp('DPC only',type)
                dpcOnly = false;
            end
        
        startNum = 1;
        Jrun = sortJrun(Jrun);
   
        %% %%%%% Process images %%%%%%
        sp = 0; skip = 0; vs = 1; spf = 0; spT = 0;

        if ZeroZ == Jrun{1,1}.Z
            % Checks if the first Z is the largest value
            mloop = startNum:min(sum(~cellfun(@isempty,Jrun),2));
        else
            % Otherwise, the loop follows this:
            mloop = min(sum(~cellfun(@isempty,Jrun),2)):-1:startNum;
        end

        Ic = cell(1,4);

            for m = mloop
    
                % Read set of 4 raw images
                for mm=1:4
    
                    fid=fopen([thisPath,'\',Jrun{mm,m}.fname]);
    
                    if (fid==-1)
                        skip=1;
                        break;
                    end
    
                    if Jrun{mm,m}.mlb
                        I0=fread(fid,[2048,2048],'uint16');
                    else
                        I0=fread(fid,[2048,2048],'uint16',0,'b');
                    end
    
                    fclose(fid);
    
                    dimsRead = size(I0);
                    % Ensure even dimensions
                    if mod(dimsRead(1),2) == 1; dimsRead(1) = dimsRead(1)-1; end
                    if mod(dimsRead(2),2) == 1; dimsRead(2) = dimsRead(2)-1; end
                    I0 = I0(1:dimsRead(1),1:dimsRead(2));
    
                    [minSize,ind] = min([dims,size(I0)]);
                    maxSize = max(size(I0));
                    if ind == 3 || ind == 4; minSize = minSize - 1; end
                    if ImOptions==1 
                        % Crop to smallest image size
                        % dims = [minSize,minSize];
                        % I0 = I0(1:dims(1),1:dims(2));                    
                        % Or pad up to largest image size
                        dims = [maxSize,maxSize];
                        R1 = ceil(maxSize/dimsRead(1));
                        R2 = ceil(maxSize/dimsRead(2));
                        I0 = repmat(I0, R1, R2); % This adds zeros to the image, so we will solve that next
                        I0(I0==0) = mean(I0(I0>0));
                        I0 = I0(1:dims(1),1:dims(2));
    
                        I0(1:5,1) = I0(1:5,2);
                        I0(6:7,2) = I0(6:7,1);
                    end
    
                    if(std(I0(:)) ~= 0)
                        Ic{mm}=I0./std(I0(:));
                    else
                        skip = 1;
                        break;
                    end
                end
    
                if skip
                    skip=0;
                    continue;
                end
    
                % This bloodProcessing function is needed to properly process
                % and display the images
                this = bloodProcessing(Ic,2,H_phase,mag,dpcOnly,type);
                if std(this{1}(:)) == 0 || isnan(std(this{1}(:)))
                    continue;
                end
    
                sp = sp + 1;
                Xpos = Jrun{1,m}.X;
                Ypos = Jrun{1,m}.Y;
    
                if sp > 1
                    if Jrun{1,m}.Z ~= Zpos && saveZs
                        if ~dpcOnly
                            toSave = matfile([OutPath,'\',fID,'_',num2str(Zpos),'.mat'],'writable',true);
                            toSave.phaseIM = phaseIM;
                            phaseIM = [];
                        end
                        toSave.positions = positions;
                        sp = 1;
                        positions = [];
                    end
                end
    
                Zpos = Jrun{1,m}.Z;
                Zdisp = (ZeroZ - Zpos);
                
                if ~dpcOnly
                    phaseim = this{5};
                    phaseim = phaseim(1:dimsRead(1),1:dimsRead(2));
                end
    
    
                if ~dpcOnly
                    phaseIM(:,:,sp) = phaseim;
                
                    thisPhase = phaseIM(51:end-50,51:end-50,sp);
                    positions(1:3,sp) = [Xpos, Ypos, Zpos];
    
                    imrgb = phaseIM(51:end-50,51:end-50,sp);
                    imrgb(imrgb<cmin) = cmin; imrgb(imrgb>cmax) = cmax;
    
                    mapped_image = (double(imrgb) - cmin) ./ (cmax - cmin);
                    mapped_image = mapped_image .* 256;
                    mapped_image = uint8(mapped_image);
                    mapped_image = ind2rgb(mapped_image,bone(256));
    
    

                        imwrite(mapped_image,...
                            [OutPath,'\images\',fID,'__',num2str(sp),'.tif'],...
                            'Compression', 'none')


    
                end            
            end
    
            if ImOptions==1 
                if ~saveZs
                    toSave = matfile([OutPath,'\',fID,'_qOBM.mat'],'writable',true);
                    toSave.phaseIM = phaseIM;
                    toSave.positions = positions;
                else
                    % If Save each Z individually is checked, then the Z
                    % position of each image is included in the file
                    toSave = matfile([OutPath,'\',fID,'_',num2str(Zpos),'.mat'],'writable',true);
                    toSave.phaseIM = phaseIM;
                    toSave.positions = positions;
                end
            end       
        end
    end
end
