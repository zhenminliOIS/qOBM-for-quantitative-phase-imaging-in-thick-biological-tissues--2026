function [mag, NA, cmin, cmax, typeTissue, fitobject, S, theta, LD, freq] = getFileInformation(fID, type, mag, defTH)

typeTissue = type;
fitobject = []; S = []; freq = []; LD = []; theta = []; cmax = []; cmin = []; NA = []; mag = []; 

% Check if File Name Includes Magnification
if contains(fID,'20x') || contains(fID,'20X')
    mag = 20;
    NA = .45;
    LD = 10;
    if strcmp('Brain',type) || strcmp('Fixed Brain',type)
        cmin = -.1;
        cmax = .35;
    else
        cmin = -.1;
        cmax = .35;
    end

elseif contains(fID,'60x') || contains(fID,'60X')
    mag = 60;
    NA = .7;
    LD = 7.25;
    if strcmp('Brain',type) || strcmp('Fixed Brain',type)
        cmin = -.4;
        cmax = 1;
    else
        cmin = -.2;
        cmax = .45;
    end

elseif (contains(fID,'40x') || contains(fID,'40X')) && (contains(fID,'NA0_95') || contains(fID,'95NA'))
    mag = 40;
    NA = .95;
    LD = 7.25;
    if strcmp('Brain',type) || strcmp('Fixed Brain',type)
        cmin = -.15;
        cmax = .25;
    else
        cmin = -.15;
        cmax = .4;
    end

elseif contains(fID,'40x') || contains(fID,'40X')
    mag = 40;
    NA = .6;
    LD = 9.4;
    if strcmp('Brain',type) || strcmp('Fixed Brain',type)
        cmin = -.4;
        cmax = 1;
    else
        cmin = -.15;
        cmax = .4;
    end
    
elseif contains(fID,'10x') || contains(fID,'10X')
    mag = 10;
    NA = .25;
    LD = 8.5;
    if strcmp('Brain',type) || strcmp('Fixed Brain',type)
        cmin = -.03;
        cmax = .12;
    else
        cmin = -.03;
        cmax = .12;
    end

else
    switch mag
        case 20
            NA = .45;
            LD = 10;
            if strcmp('Brain',type) || strcmp('Fixed Brain',type)
                cmin = -.1;
                cmax = .35;
            else
                cmin = -.1;
                cmax = .35;
            end
        case 60
            NA = .7;
            LD = 7.25;
            if strcmp('Brain',type) || strcmp('Fixed Brain',type)
                cmin = -.4;
                cmax = 1;
            else
                cmin = -.2;
                cmax = .45;
            end
        case 40
            NA = .6;
            LD = 9.4;
            if strcmp('Brain',type) || strcmp('Fixed Brain',type)
                cmin = -.4;
                cmax = 1;
            else
                cmin = -.15;
                cmax = .4;
            end
        case 10
            NA = .25;
            LD = 8.5;
            if strcmp('Brain',type) || strcmp('Fixed Brain',type)
                cmin = -.1;
                cmax = .35;
            else
                cmin = -.1;
                cmax = .35;
            end
        end
end


% Check if File Name Includes Fiber Angle
if contains(fID,'65deg') || contains(fID,'65d')
    theta = 65;
    LD = 5; % Horizontal distance between fibers and objective
    if strcmp('Brain',type)
        load('MC\Source_BRN_65_720.mat','S','fitobject');
    else
        disp('WARNING - No transfer function saved for this combination')
        typeTissue = 'Paper';
    end

elseif contains(fID,'55deg') || contains(fID,'55d')
    theta = 55;
    LD = 7.25; % Horizontal distance between fibers and objective
    if strcmp('Brain',type)
        load('MC\Source_BRN_55deg_60X_LD5_Hei0_720nm.mat','out');
        S=out.S;
        fitobject=out.fitobject;
    elseif strcmp('Fixed Brain',type)
        load('MC\Source_BRNFIX_55_farRed.mat','fitobject','S'); 
    else
        disp('WARNING - No transfer function saved for this combination')
        typeTissue = 'Paper';
    end

elseif contains(fID,'45deg') || contains(fID,'45d')
    theta = 45;
    %LD = 7; % Horizontal distance between fibers and objective
    if strcmp('Brain',type)
        load('MC\Source_BRN_45_farRed.mat','fitobject','S');
    elseif strcmp('Fixed Brain',type)
        load('MC\Source_BRNFIX_45_farRed_2.mat','fitobject','S');
    else
        disp('WARNING - No transfer function saved for this combination')
        typeTissue = 'Paper';
    end

else
    theta = defTH;
    if defTH == 45
        %LD = 7; % Horizontal distance between fibers and objective
        if strcmp('Brain',type)
            load('MC\Source_BRN_45_farRed.mat','fitobject','S');
        elseif strcmp('Fixed Brain',type)
            load('MC\Source_BRNFIX_45_farRed_2.mat','fitobject','S');
        end
    elseif defTH == 65
        LD = 5; % Horizontal distance between fibers and objective
        if strcmp('Brain',type)
            load('MC\Source_BRN_65_720.mat','S','fitobject');
        else
            disp('WARNING - No transfer function saved for this combination')
            typeTissue = 'Paper';
        end
    end
end

