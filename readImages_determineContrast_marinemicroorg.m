clear
clc

% Pixel per µm in the microscope images 
pixperDist = 10.96; % per µm
distperPix = 1/pixperDist; % distance per pixel


% Loading images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filefolder = '../Data/';
filename = 'helix2_DF_1417710us_gain20.bmp';
filenameBF = 'helix2_BF_3006us_gain0.bmp';

% read raw image files for SLED and BF
imageRaw = imread(char(strcat(filefolder, filename)));
imageRawBF = imread(char(strcat(filefolder, filenameBF)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Cropping images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define center of area to crop in the images
CropCenter = [910; 1800];
CropCenterBF = [910; 1800];

% Crop dimensions 
ImWidth = 650;
ImHeight = 500;
Dimensions = [ImHeight; ImWidth];


% crop image to relevant areas 
imageCrop = imageRaw(CropCenter(1)-floor(Dimensions(1)/2):CropCenter(1)+floor(Dimensions(1)/2),...
                      CropCenter(2)-floor(Dimensions(2)/2):CropCenter(2)+floor(Dimensions(2)/2),...
                        :);
                    
imageCropBF = imageRawBF(CropCenterBF(1)-floor(Dimensions(1)/2):CropCenterBF(1)+floor(Dimensions(1)/2),...
                      CropCenterBF(2)-floor(Dimensions(2)/2):CropCenterBF(2)+floor(Dimensions(2)/2),...
                        :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Background removal and salt/pepper filter (noise on camera) %%%%%%%%%%%%%
y1 = 2;
y2 = 4;
x1 = 2;
x2 = 4; 

bkgd = floor(sum(sum(imageCrop(y1:y2,x1:x2,1)))/(x2-x1)/(y2-y1));
imageCropBkgd = imageCrop - bkgd;

image = imadjust(imageCropBkgd, [0 0 0; 0.7 1 1], []);
imageRed = medfilt2(image(:,:,1));

image(:,:,1) = imageRed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Cut line across micro organism %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = 295;
x2 = 396;
y1 = 188;
y2 = 288;

%along micro organism
% x1 = 272;
% x2 = 377;
% y1 = 269;
% y2 = 194;

buffer = 3;
hbar = 3;
barlength = 5*pixperDist;

figure()
  imshow(imageCrop)
  line([x1 x2], [y1 y2], 'LineStyle', ':', 'Linewidth', 1.5, 'Color', [1 1 1]);
  rectangle('Position', [ImWidth-buffer-barlength ImHeight-buffer-hbar barlength hbar], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
  
figure()
  imshow(imageCropBF(:,:,2))    
  line([x1 x2], [y1 y2], 'LineStyle', ':', 'Linewidth', 1.5, 'Color', [0 0 0]); 
  rectangle('Position', [ImWidth-buffer-barlength ImHeight-buffer-hbar barlength hbar], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


% Extract intensity profile %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cx,cy,DFprofile] = improfile(imageCrop,[x1 x2],[y1 y2]);  

dx = abs(x2-x1);
dy = abs(y2-y1);
xProfile = distperPix*sqrt(dx^2+dy^2)/size(cx,1)*(0:1:size(cx,1)-1);

n = 2; 
DFprofiles = zeros(size(cx,1),2*n+1);
DFprofilesBF = zeros(size(cx,1),2*n+1);

for ii=-n:1:n
    profile = improfile(imageCrop,[x1-ii x2-ii],[y1-ii y2-ii]); 
    profileBF = improfile(imageCropBF,[x1-ii x2-ii],[y1-ii y2-ii]); 
    DFprofiles(:,ii+n+1) = squeeze(profile(:,1,1));
    DFprofilesBF(:,ii+n+1) = squeeze(profileBF(:,1,2)); %green channel
end

DFprofileAv = mean(DFprofiles,2);
DFprofileAv = DFprofileAv/max(DFprofileAv);
DFprofileAvBF = mean(DFprofilesBF,2);
DFprofileAvBF = DFprofileAvBF/max(DFprofileAvBF);

% fraction of data points above high threshold and below low threshold for
% contrast determination 
fraction = 0.15;

% contrast from intensity profile 
[CDF, dCDF, nNumforAv, minMaxGrDF, maxMinGrDF] = findContrast(DFprofileAv,fraction)
[CBF, dCBF, nNumforAvBF, minMaxGrBF, maxMinGrBF] = findContrast(DFprofileAvBF,fraction)     


% plot intensity profile and thresholds used for contrast determination %%%
 
limoffset = 0.05;

capSize = 9;
numSize = 9;
        
figure()
    hAxis(1) = subplot(2,1,1)
        hold on 
        plot(xProfile, DFprofileAvBF, 'Linewidth', 1, 'Color', [0 0 0])
        line([min(xProfile) max(xProfile)], [maxMinGrBF maxMinGrBF], 'Linestyle', ':', 'Linewidth', 0.5, 'Color', [0 0 0]);
        line([min(xProfile) max(xProfile)], [minMaxGrBF minMaxGrBF], 'Linestyle', ':', 'Linewidth', 0.5, 'Color', [0 0 0]);
        
        set(gca, 'Fontsize', numSize)
        xlim([min(xProfile) max(xProfile)])
        ylim([min(DFprofileAvBF)-limoffset max(DFprofileAvBF)+limoffset])
        xticklabels([])
        box on
        
        ylBF = ylim
        
        
    hAxis(2) = subplot(2,1,2)
        plot(xProfile, DFprofileAv, 'Linewidth', 1, 'Color', [1 0 0])
        line([min(xProfile) max(xProfile)], [maxMinGrDF maxMinGrDF], 'Linestyle', ':', 'Linewidth', 0.5, 'Color', [0 0 0]);
        line([min(xProfile) max(xProfile)], [minMaxGrDF minMaxGrDF], 'Linestyle', ':', 'Linewidth', 0.5, 'Color', [0 0 0]);        
        
        %ylim([min(DFprofileAv)-limoffset min(DFprofileAv)-limoffset+limdelta])
        xlim([min(xProfile) max(xProfile)])
        ylim([min(DFprofileAv)-limoffset max(DFprofileAv)+limoffset])
        ylDF = ylim
        
        
%         NumTicks = 7;
%         L = get(gca,'XLim');
%         L(2) = floor(L(2)/NumTicks)*NumTicks;
        set(gca,'XTick',[0 2 4 6 8 10 12 14])
        
        set(gca, 'Fontsize', numSize)
        set(gcf, 'Color', [1 1 1])
        box on       
       
        xlabel('distance (\mum)', 'Fontsize', capSize)
    
        sizeYDF = ylDF(2) - ylDF(1);
        sizeYBF = ylBF(2) - ylBF(1);
        
        pad = 0.05;
        bottom = 0.2;
        
        sizeTotal = sizeYDF + sizeYBF;
        sizeTotalPad = (1+bottom+2*pad)*sizeTotal;
        
        pos = get( hAxis(1), 'Position' )
        pos(2) = (sizeYDF + (pad+bottom)*sizeTotal)/ sizeTotalPad ;                         % Shift down.
        pos(4) = sizeYBF/sizeTotalPad ;                        % Increase height.
        set( hAxis(1), 'Position', pos ) 
        
        pos = get( hAxis(2), 'Position' )
        pos(2) = (bottom*sizeTotal) / sizeTotalPad ;                         % Shift down.
        pos(4) = sizeYDF/sizeTotalPad ;                        % Increase height.
        set( hAxis(2), 'Position', pos ) 
        
             
        p1=get(hAxis(1),'position');
        p2=get(hAxis(2),'position');

        height= p1(2)+ p1(4)-p2(2);
        h3=axes('position',[0.65*p2(1) p2(2) p2(3) height],'visible','off');
       
        h_label=ylabel('normalized intensity','visible','on', 'Fontsize', capSize);

        set(gcf, 'Color', [1 1 1]) 
        %set(gcf, 'Menubar', 'none', 'Toolbar', 'none') 
        
        set(gcf, 'Units', 'points', 'Position', [50 50 140 95])
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
