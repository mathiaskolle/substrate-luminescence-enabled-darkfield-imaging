
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Calculations for partial coherent imaging with SEDFI surfaces and     %
%   bright field illumination                                             %
%   - Phase coherence factor                                              %
%   - Amplitude transfer fct. and amplitude spread fct. of 4f-system      %
%   - Phase distribution of colloidal object                              %
%                                                                         %
%     Mathias Kolle - Oct. 31, 2019                                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% Microscope properties and other boundary conditions  
refrInd = 1.33;  % refractive index of surrounding medium (imaging medium, here water)
ObjNA = 1; % NA of objective 
lambda0 = 630E-9; % wavelength of illumination in meters

fcObj = ObjNA / lambda0; % cut off in the ATF of the objective
fmax = multiplier * fcObj; % maximum frequency for which to evaluate ATF (much larger than the cut-off to get finer resolution when doing Fourier transform) 


% Light angles   
thMin = 64;
thMax = 85;

%plot parameters 
sizeFont = 24; % font size for plots
ASFplotCrop = 0.15; % factor that defines for what spatial extend the amplitude spread function data is shown in plots
TwoDimSize = 0.25; %size of data images for ATF and ASF 
ATFplotCrop = 0.1; % ratio that defines for what spatial frequency range the ATF data is shown in plots
multiplier = 20; % used to define frequency range for ATF data


dx = 1/2/fmax; % stepwidth in real space where the ASF lives 
nPoints = 1000; % number of points in the frequency range from -fmax to fmax
deltaf = 1/nPoints/dx; % frequency resolution


% variables used to define scale bar in ATF plots 
ATFbuffer = 0.3;  
ATFbheight = 0.2;
ATFbarlength = 1;

ASFbuffer = 0.7E-7;  
ASFbheight = 1E-7;
ASFbarlength = 1E-6;


%%%%%%%%%%%%%%%%%% Amplitude transfer functions %%%%%%%%%%%%%%%%%%%%%%%%%%%

farray = -fmax:deltaf:fmax;
nSize = size(farray,2);

fx = farray;
fy = farray';
fr = sqrt(fx.^2 + fy.^2); % since the aperture has circular symetry, work with radial coordinates 

% ATF plot limits
limfx = ATFplotCrop*max(fx)*1E-6;

% ATF of Objective
ATFObj = 1*(fr<fcObj); 

% plot ATF of Objective
figure('Units', 'normalized', 'Position', [0.1 0.1 TwoDimSize TwoDimSize])
pcolor(fx*1E-6,fy*1E-6, ATFObj)
xlim([-limfx limfx])
ylim([-limfx limfx])
shading flat
rectangle('Position', [limfx-ATFbuffer-ATFbarlength, -limfx+ATFbuffer, ATFbarlength, ATFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
axis square
xticks([])
yticks([])
set(gca, 'Fontsize', sizeFont)
c1 = colorbar
c1.Label.String = 'Amplitude';
set(c1,'YTick',-1:0.2:1)
colormap gray
set(gcf, 'Color', [1 1 1])


%onset and cut-off frequencies of light source - here simplified to a
% binary annulus (the more realistic ATF for the source is formed later below 
fclowLight = refrInd * sind(thMin)/lambda0;
fchighLight = refrInd * sind(thMax)/lambda0;
    
% ATF of light source 
ATFLight = 1*(fr<=fchighLight & fr>fclowLight); 

% plot ATF of light source
figure('Units', 'normalized', 'Position', [0.1 0.1 TwoDimSize TwoDimSize])
pcolor(fx*1E-6,fy*1E-6, ATFLight)
xlim([-limfx limfx])
ylim([-limfx limfx])
shading flat
rectangle('Position', [limfx-ATFbuffer-ATFbarlength, -limfx+ATFbuffer, ATFbarlength, ATFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
axis square
xticks([])
yticks([])
set(gca, 'Fontsize', sizeFont)
c1 = colorbar
c1.Label.String = 'Amplitude';
set(c1,'YTick',-1:0.2:1)
colormap gray
set(gcf, 'Color', [1 1 1])


%%%%%%%%%%%%%%%%%% Amplitude spread functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% space co-ordinates for Amplitude Spread Functions (ASF)
xASFObj = -1/2/deltaf:dx:1/2/deltaf;
yASFObj = xASFObj';

% Amplitude spread functions for Objective and light source from the
% Fourier transforms of their ATFs 
%ASFObj = fftshift(fft2(ATFObj, nSize, nSize))/sqrt(length(fx)*length(fy)); 
ASFObj = fftshift(fft2(ATFObj))/sqrt(length(fx)*length(fy)); 
ASFLight = fftshift(fft2(ATFLight))/sqrt(length(fx)*length(fy)); 

% plot limit for ASFs 
axlim = ASFplotCrop*max(xASFObj);

% Making colormap for ASFs
colormap gray
cmap = colormap;    
cmapflip = flipud(cmap);
cmapflip(:,1) = 1;
map = [cmap; cmapflip];
    
    
% ASF can also be directly obtained with the theoretical expression
% involvong Bessel fct.; if done right, it matches the ASFs obtained above
% through Fourier transformation of the ATFs


% ASF for objective
besselarg = 2*pi/lambda0 * ObjNA * sqrt(xASFObj.^2 + yASFObj.^2);
ASFObjTheo = 2 * besselj(1, besselarg) ./ (besselarg); % coherent     
ASFObjTheo(isnan(ASFObjTheo)) = 1;
%ASFObjTheo = pi*fcObj^2/(4*fmax^2)*nPoints*ASFObjTheo; % normalize to match with Fourier transform of ATF
ASFObjTheo = sum(sum(abs(ATFObj).^2))/nPoints*ASFObjTheo;



figure('Units', 'normalized', 'Position', [0.1 0.1 TwoDimSize TwoDimSize])
    pcolor(xASFObj, yASFObj, ASFObjTheo)    
    rectangle('Position', [axlim-ASFbuffer-ASFbarlength, -axlim+ASFbuffer, ASFbarlength, ASFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
    xlim([-axlim axlim])
    ylim([-axlim axlim])
    shading flat    
    axis square
    xticks([])
    yticks([])
    set(gca, 'Fontsize', sizeFont)
    c1 = colorbar;
    c1.Label.String = 'Amplitude';  
    set(c1,'YTick',-1:0.2:1)
    caxis([-0.5 1])
    colormap(map)
    set(gcf, 'Color', [1 1 1])

    
checkersY = -floor(size(ASFObj,1)/2):1:floor(size(ASFObj,1)/2);
checkersX = checkersY';
checkers = (2*mod(checkersX,2)-1)*(2*mod(checkersY,2)-1);
ASFObjCor = checkers.*sign(real(ASFObj)).*abs(ASFObj);
ASFLightCor = checkers.*sign(real(ASFLight)).*abs(ASFLight);  
    
    
% Theoretical ASF for light source = normalised difference between two ASFs fof maximum light
% incidence angle and minimum light incidence angle 
ASFthMin = 2 * besselj(1, 2*pi/lambda0 * refrInd * sind(thMin) * sqrt(xASFObj.^2 + yASFObj.^2)) ./ (2*pi/lambda0 * refrInd * sind(thMin) * sqrt(xASFObj.^2 + yASFObj.^2)) ;
ASFthMax = 2 * besselj(1, 2*pi/lambda0 * refrInd * sind(thMax) * sqrt(xASFObj.^2 + yASFObj.^2)) ./ (2*pi/lambda0 * refrInd * sind(thMax) * sqrt(xASFObj.^2 + yASFObj.^2)) ;
ASFthMin(isnan(ASFthMin)) = 1;    
ASFthMax(isnan(ASFthMax)) = 1;
ASFLightTheo = 1 / ( sind(thMax)^2 - sind(thMin)^2 ) * ( sind(thMax)^2 * ASFthMax - sind(thMin)^2 * ASFthMin );
ASFLightTheo = sum(sum(abs(ATFLight).^2))/nPoints*ASFLightTheo;


figure('Units', 'normalized', 'Position', [0.1 0.1 0.4 0.6])
    %hAxis(1) = subplot(2,1,1)
        hold on    
        plot(xASFObj*1E6, ASFObjTheo(:,round(size(ASFObjTheo,2)/2)), 'Linewidth', 1.5, 'Color', [0 0 0]) 
        plot(xASFObj*1E6, ASFObjCor(:,round(size(ASFObjCor,2)/2)), 'Linewidth', 1.5, 'Color', [1 0 0])
       
        plot(xASFObj*1E6, ASFLightCor(:,round(size(ASFLightCor,2)/2)), 'Linewidth', 1.5, 'Color', [0 0 1]) 
        plot(xASFObj*1E6, ASFLightTheo(:,round(size(ASFLightTheo,2)/2)), 'Linewidth', 1.5, 'Color', [0 1 0]) 
        
        line([-axlim*1E6 axlim*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        xlim([-axlim*1E6 axlim*1E6])
        ylim([min(min(ASFObjCor))-0.1*max(max(ASFObjCor)) 1.1*max(max(ASFObjCor))])
        set(gca, 'Fontsize', sizeFont)
        xticks([])
        box on;   


figure()
    pcolor(xASFObj, yASFObj, ASFLightTheo)
    rectangle('Position', [axlim-ASFbuffer-ASFbarlength, -axlim+ASFbuffer, ASFbarlength, ASFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
    xlim([-axlim axlim])
    ylim([-axlim axlim])
    shading flat    
    axis square
    xticks([])
    yticks([])
    set(gca, 'Fontsize', sizeFont)
    c1 = colorbar
    c1.Label.String = 'Amplitude';   
    colormap(map)
    set(gcf, 'Color', [1 1 1])
   


% Obtaining ATF of light source from theoretical ASF
ATFLightTheo = fftshift(fft2(ASFLightTheo, nSize, nSize)); 

figure()
    pcolor(fx*1E-6,fy*1E-6, abs(ATFLightTheo)/max(max(abs(ATFLightTheo))))
    xlim([-limfx limfx])
    ylim([-limfx limfx])
    shading flat
    rectangle('Position', [limfx-ATFbuffer-ATFbarlength, -limfx+ATFbuffer, ATFbarlength, ATFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
    axis square
    xticks([])
    yticks([])
    set(gca, 'Fontsize', sizeFont)
    c1 = colorbar
    c1.Label.String = 'Amplitude';
    colormap gray
    set(gcf, 'Color', [1 1 1])


%%%%%%% Calculate ATF of light with accounting for Bragg reflector %%%%%%%
% use calculated intensity vs. angle distribution from Bragg reflector
% Gold reflection
Data = dlmread('Au_refr_ind_good.csv');
A_lambda = Data(:,1)*1000;
ng = Data(:,2);
kg = Data(:,3);

% Quantum dot emission 
Data = dlmread('ic_i_39_emission.csv');
qd_em_lambda = Data(:,1);
qd_em = Data(:,2); 

% Quantum dot absorption 
Data = dlmread('ic_i_39_abs.csv');
qd_abs_lambda = Data(:,1);
qd_abs = Data(:,2);

% Interpolation for getting finer sampling
sample_num = 1000;
lambda = (min(A_lambda):(max(A_lambda)-min(A_lambda))/(sample_num-1):max(A_lambda))';
ng_fine = interp1(A_lambda,ng,lambda, 'pchip', 0);
kg_fine = interp1(A_lambda,kg,lambda, 'pchip', 0);

qd_em_fine = interp1(qd_em_lambda,qd_em,lambda, 'pchip', 0);
qd_abs_fine = interp1(qd_abs_lambda,qd_abs,lambda, 'pchip', 0);



lambda0tc = qd_em_lambda(find(qd_em == max(qd_em)))
min(qd_em)
hFWHM = (max(qd_em) - min(qd_em))/2

lambdaFWHM = find((qd_em - hFWHM).^2 < 1E-3)
lambda(lambdaFWHM(1))
lambda(lambdaFWHM(end))
dlambdatc = lambda(lambdaFWHM(end)) - lambda(lambdaFWHM(1))

c = 3E8;
tc = (lambda0tc*1E-9)^2/(c * dlambdatc*1E-9)

d = 50E-6;
dn = 0.02; 
tdelay = dn * d / c





%hFWHM = 


tc = lambda0tc^2


sizeFont = 24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialisation of the relevant parameters
I = complex(0,1);


%refractive index layer 1 and layer 2
n1 = 2.2;
n2 = 1.49;

% center wavelength of Bragg reflector 1 and 2
l_center1 = 630;%630;%585;%610;
l_center2 = 550;

% thicknesses in Bragg reflectors
d1 = l_center1/4/n1; 
d2 = l_center1/4/n2;

d3 = l_center2/4/n1;
d4 = l_center2/4/n2;

%Layer numbers in the two Bragg reflectors
nl1 = 13;
nl2 = 0;

% index of medium outside of sample
nOutside = refrInd;%[1; 1.33; 1.518];
% index of refraction of upper medium (where coming from - here PLMA or PMMA)
n_up = 1.49;

ai = (0:0.5:asind(refrInd/n_up))';  %60x water

% preallocation of arrays for layer thicknesses and refractive indices 
n = zeros(nl1+nl2,1);
d = zeros(nl1+nl2,1);

% construction of stacks (thicknesses and refractive indices
for ii=1:1:nl1+nl2
   n(ii) = (ii<=nl1)*(mod(ii,2)*n1+mod(ii-1,2)*n2) + (ii > nl1)*(mod(ii-nl1,2)*n1+mod(ii-nl1-1,2)*n2);
   d(ii) = (ii<=nl1)*(mod(ii,2)*d1+mod(ii-1,2)*d2) + (ii > nl1)*(mod(ii-nl1,2)*d3+mod(ii-nl1-1,2)*d4);
end

% Data matrices 
Rp_angles = zeros(size(lambda,1), size(ai,1), size(nOutside,1));
Rs_angles = zeros(size(lambda,1), size(ai,1), size(nOutside,1));
for jj = 1:1:size(nOutside,1)
    n_low = nOutside(jj)+zeros(size(ng_fine,1),1);
    for ii = 1:1:size(ai,1)
        x0 = [ai(ii); n_up];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate Reflectivities for all angles  (d, n, x0, n_low, lambda);
        [Rp,Rs] = multilayer_refl_v3(d, n, x0, n_low, lambda);
        Rp_angles(:, ii, jj) = Rp;
        Rs_angles(:, ii, jj) = Rs;       
    end 
end

Rtot_angles = (Rp_angles + Rs_angles)/2; % unpolarized light 


% full signal taking into account QD emission
signal = (1-Rtot_angles).*repmat(qd_em_fine, 1, size(Rtot_angles,2), size(Rtot_angles,3));

% emission angles 
ai_out = asind(n_up/nOutside.*sind(ai));

% find matching wavelength (630nm in Bragg reflector data)
index = find(abs(lambda-lambda0*1E9)<1);
spfreq = sind(ai_out)*nOutside/(lambda(index(1))*1E-9); %deleted the x2 here 07/16/19
signalLambda0 = signal(index(1),:,1);% - signal(index(1),1,1);
signalLambda0(end) = 0;

    
% build 2D ATF of light source 
signalLambda0twoWings = [fliplr(signalLambda0(2:end)), signalLambda0];  
spfreqtwoWings = [-flipud(spfreq(2:end)); spfreq];

% figure() 
%     plot(spfreqtwoWings*1E-6, signalLambda0twoWings);
%     hold on
%     plot(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)));
%     xlim([-3 3])

% evaluate ATF of light source at frequencies for which Objective ATF was calculated    
signalLambda0twoWingsfx = interp1(spfreqtwoWings,signalLambda0twoWings,fx, 'pchip', 0);    
figure() 
    hold on
    area(fx*1E-6, signalLambda0twoWingsfx/max(signalLambda0twoWingsfx), 'Linewidth', 1.5, 'Facecolor', [1 0.8 0.8], 'Edgecolor', 'none');
    area(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Facecolor', [0.8 0.8 0.8], 'Edgecolor', 'none', 'FaceAlpha', 0.5);
    plot(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Color', [0 0 0]);
    plot(fx*1E-6, signalLambda0twoWingsfx/max(signalLambda0twoWingsfx), 'Linewidth', 1.5, 'Color', [1 0 0]);
    
    xlim([-limfx limfx])
    ylim([-0.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx))) 1.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx)))])
    set(gca, 'Fontsize', sizeFont)
    xlabel('cycles 1/\mum');
    ylabel('Magnitude')
    set(gcf,'Color', [1 1 1]);
    box on;
   

% build 2D ATF 
fx_mat = zeros(size(fx,2),size(fx,2));
Illum_mat = zeros(size(fx,2),size(fx,2));
for ii = 1:1:size(fx_mat,1)
    for jj = 1:1:size(fx_mat,2)
       fx_mat(ii,jj) = sqrt(fx(ii)^2 + fx(jj)^2);
       [M,Ind] = min(abs(fx - fx_mat(ii,jj)));
        Illum_mat(ii,jj) = signalLambda0twoWingsfx(1,Ind,1);
    end
end    
    

% plot 2D ATF of light source 
figure()
    pcolor(fx*1E-6,fy*1E-6,  Illum_mat)
    xlim([-limfx limfx])
    ylim([-limfx limfx])
    shading flat
    rectangle('Position', [limfx-ATFbuffer-ATFbarlength, -limfx+ATFbuffer, ATFbarlength, ATFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
    axis square
    xticks([])
    yticks([])
    set(gca, 'Fontsize', sizeFont)
    c1 = colorbar
    c1.Label.String = 'Amplitude';
    colormap gray
    set(gcf, 'Color', [1 1 1])
        
        

%%%% Amplitude spread function for angle spectrum of SEDFI source %%%%%%%%%
ASFrealLight = fftshift(fft2(Illum_mat, nSize, nSize))/sqrt(length(fx)*length(fy)); ; 

% need to take care of the 1-pixel phase undulation in the Fourier transform of of the ATF in Illum_mat ! 
checkersY = -floor(size(ASFrealLight,1)/2):1:floor(size(ASFrealLight,1)/2);
checkersX = checkersY';
checkers = (2*mod(checkersX,2)-1)*(2*mod(checkersY,2)-1);
ASFrealLightCor = checkers.*sign(real(ASFrealLight)).*abs(ASFrealLight);

% plot ASF of SEDFI light source
figure()
    pcolor(xASFObj, yASFObj, ASFrealLightCor)
    rectangle('Position', [axlim-ASFbuffer-ASFbarlength, -axlim+ASFbuffer, ASFbarlength, ASFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
    xlim([-axlim axlim])
    ylim([-axlim axlim])
    shading flat    
    axis square
    xticks([])
    yticks([])
    set(gca, 'Fontsize', sizeFont)
    c1 = colorbar
    c1.Label.String = 'Amplitude'; 
    colormap(map)
    set(gcf, 'Color', [1 1 1])


figure('Units', 'normalized', 'Position', [0.1 0.1 0.4 0.6])
    hAxis(1) = subplot(2,1,1)
        plot(fx*1E-6, signalLambda0twoWingsfx, 'Linewidth', 1.5, 'Color', [1 0 0]);
        xlim([-limfx limfx])
        ylim([-0.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx))) 1.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx)))])
        set(gca, 'Fontsize', sizeFont) 
        xticks([])
        box on;
        
    hAxis(2) = subplot(2,1,2)
        plot(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Color', [0 0 0]);
        xlim([-limfx limfx])
        ylim([-0.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx))) 1.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx)))])
        set(gca, 'Fontsize', sizeFont)        
        xlabel('cycles 1/\mum');
        box on;
    
        pad = 0.1;
        bottom = 0.4;
        
        sizeTotalPad = (2+bottom+2*pad);

        pos = get( hAxis(1), 'Position' )
        pos(2) = (1 + pad + bottom)/ sizeTotalPad ;                         % Shift down.
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(1), 'Position', pos ) 

        pos = get( hAxis(2), 'Position' )
        pos(2) = bottom / sizeTotalPad ;                         % Shift down.
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(2), 'Position', pos )

        
        p1=get(hAxis(1),'position');
        p2=get(hAxis(2),'position');

        height=p1(2)+p1(4)-p2(2);
        h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
        set(gca, 'Fontsize', sizeFont)
        h_label=ylabel('magnitude','visible','on');

        set(gcf, 'Color', [1 1 1])
    
limYmax = 1.1 * max( max(max(ASFrealLightCor)), max(max(ASFObjTheo)) );
limYmin = min(min(min(ASFrealLightCor)), min(min(ASFObjTheo))) - 0.1/1.1*limYmax;     

figure('Units', 'normalized', 'Position', [0.1 0.1 0.4 0.6])
    hAxis(1) = subplot(2,1,1)
        %plot(xASFObj*1E6, ASFObjTheo(:,round(size(ASFObjTheo,2)/2)), 'Linewidth', 1.5, 'Color', [0 0 0]) 
        %plot(xASFObj*1E6, ASFLightTheo(:,round(size(ASFLightTheo,2)/2)), 'Linewidth', 1.5, 'Color', [0 0 1]) 
        plot(xASFObj*1E6, ASFrealLightCor(:,round(size(ASFrealLightCor,2)/2)), 'Linewidth', 1.5, 'Color', [1 0 0])
        hold on
        plot(xASFObj*1E6, ASFLightTheo(:,round(size(ASFLightTheo,2)/2)), 'Linewidth', 1.5, 'Color', [0 1 0])
        plot(xASFObj*1E6, ASFLightCor(:,round(size(ASFLightCor,2)/2)), 'Linewidth', 1.5, 'Color', [0 0 1])
        
        line([-axlim*1E6 axlim*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        xlim([-axlim*1E6 axlim*1E6])
        ylim([limYmin limYmax])
        set(gca, 'Fontsize', sizeFont)
        xticks([])
        box on;
        
    hAxis(2) = subplot(2,1,2)
        plot(xASFObj*1E6, ASFObjTheo(:,round(size(ASFObjTheo,2)/2)), 'Linewidth', 1.5, 'Color', [0 0 0]) 
        %plot(xASFObj*1E6, ASFLightTheo(:,round(size(ASFLightTheo,2)/2)), 'Linewidth', 1.5, 'Color', [0 0 1]) 
        %plot(xASFObj*1E6, ASFrealLightCor(:,round(size(ASFrealLightCor,2)/2)), 'Linewidth', 1.5, 'Color', [1 0 0])
        line([-axlim*1E6 axlim*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        xlim([-axlim*1E6 axlim*1E6])
        ylim([limYmin limYmax])
        set(gca, 'Fontsize', sizeFont)
        xlabel('distance (\mum)');
        box on;
    
        pad = 0.1;
        bottom = 0.4;
        
        sizeTotalPad = (2+bottom+2*pad);

        pos = get( hAxis(1), 'Position' )
        pos(2) = (1 + pad + bottom)/ sizeTotalPad ;                         % Shift down.
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(1), 'Position', pos ) 

        pos = get( hAxis(2), 'Position' )
        pos(2) = bottom / sizeTotalPad ;                         % Shift down.
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(2), 'Position', pos )

        
        p1=get(hAxis(1),'position');
        p2=get(hAxis(2),'position');

        height=p1(2)+p1(4)-p2(2);
        h3=axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
        set(gca, 'Fontsize', sizeFont)
        h_label=ylabel('normalized amplitude','visible','on');
        %ylabel('Normalized intensity');

        set(gcf, 'Color', [1 1 1])
    
% plotting ASF for paper figure

colormap gray
cmap = colormap;    
cmaptop = cmap;
cmaptop(:,1) = 1;

cmapbottom = cmap;
cmapbottom(:,2) = 0;
cmapbottom(:,3) = 0;
map = [cmapbottom; cmaptop];

ASFrealLightCorNorm = ASFrealLightCor/max(max(ASFrealLightCor));
ASFObjTheoNorm = ASFObjTheo/max(max(ASFObjTheo));

limYmax = 1.1 * max( max(max(ASFrealLightCorNorm)), max(max(ASFObjTheoNorm)) );
limYmin = min(min(min(ASFrealLightCorNorm)), min(min(ASFObjTheoNorm))) - 0.1/1.1*limYmax;     

limYmaxMap = max( max(max(ASFrealLightCorNorm)), max(max(ASFObjTheoNorm)) );
limYminMap = min( min(min(ASFrealLightCorNorm)), min(min(ASFObjTheoNorm)) );
%%
figure('Units', 'normalized', 'Position', [0.1 0.1 0.4 0.6])
    hAxis(1) = subplot(2,2,1)
        pcolor(xASFObj, yASFObj, ASFrealLightCorNorm)
        rectangle('Position', [axlim-ASFbuffer-ASFbarlength, -axlim+ASFbuffer, ASFbarlength, ASFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
        xlim([-axlim axlim])
        ylim([-axlim axlim])
        shading flat    
        axis square
        xticks([])
        yticks([])
        set(gca, 'Fontsize', sizeFont)
        %c1 = colorbar
        %c1.Label.String = 'Amplitude'; 
        colormap(map)
        caxis([limYminMap limYmaxMap])
        set(gcf, 'Color', [1 1 1])
    
    hAxis(2) = subplot(2,2,2)
        pcolor(xASFObj, yASFObj, ASFObjTheoNorm)
        rectangle('Position', [axlim-ASFbuffer-ASFbarlength, -axlim+ASFbuffer, ASFbarlength, ASFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
        xlim([-axlim axlim])
        ylim([-axlim axlim])
        shading flat    
        axis square
        xticks([])
        yticks([])
        
        c1 = colorbar
        c1.Label.String = ({'normalized'; 'amplitude'})
        set(gca, 'Fontsize', sizeFont)
        colormap (hAxis(2), gray)
        caxis([limYminMap limYmaxMap])
        set(gcf, 'Color', [1 1 1])
    
    hAxis(3) = subplot(2,2,[3 4])
        hold on
        area(xASFObj*1E6, ASFObjTheoNorm(:,round(size(ASFObjTheo,2)/2)), 'Linewidth', 1.5, 'Facecolor', [0.8 0.8 0.8], 'Edgecolor', 'none', 'FaceAlpha', 0.5)       
        plot(xASFObj*1E6, ASFObjTheoNorm(:,round(size(ASFObjTheo,2)/2)), 'Linewidth', 1.5, 'Color', [0 0 0]) 
        area(xASFObj*1E6, ASFrealLightCorNorm(:,round(size(ASFrealLightCor,2)/2)), 'Linewidth', 1.5, 'Facecolor', [1 0.8 0.8], 'Edgecolor', 'none', 'FaceAlpha', 0.5)
        plot(xASFObj*1E6, ASFrealLightCorNorm(:,round(size(ASFrealLightCor,2)/2)), 'Linewidth', 1.5, 'Color', [1 0 0])
         
        line([-axlim*1E6 axlim*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        xlim([-axlim*1E6 axlim*1E6])
        ylim([limYmin limYmax])
        xticks([-2 -1.5 -1, -0.5, 0, 0.5, 1 1.5 2])
        set(gca, 'Fontsize', sizeFont)
        xlabel({'x_2 - x_1 (\mum) @ y_2 - y_1 = 0'});
        ylabel({'normalized'; 'amplitude'},'visible','on');
        
        box on;
            
        pad = 0.1;
        padHor = 0.06;
        bottom = 0.4;
        
        sizeTotalPad = (2+bottom+2*pad);

        pos = get( hAxis(3), 'Position' )
        pos(2) = bottom/ sizeTotalPad ;                         % Shift down.
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(3), 'Position', pos ) 

        pos = get( hAxis(1), 'Position' )
        pos(1) = padHor;
        pos(2) = (1 + pad + bottom)/ sizeTotalPad ; 
        pos(3) = 0.45; % Shift down.
        pos(4) = 0.85/sizeTotalPad ;                        % Increase height.
        set( hAxis(1), 'Position', pos )

        pos = get( hAxis(2), 'Position' )
        pos(1) = 0.5-2*padHor;
        pos(2) = (1 + pad + bottom)/ sizeTotalPad ; 
        pos(3) = 0.45; % Shift down.
        pos(4) = 0.85/sizeTotalPad ;                        % Increase height.
        set( hAxis(2), 'Position', pos )  

        set(gcf, 'Color', [1 1 1])
    

% plotting ATF for paper figure   

limYmax = 1.1*max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx));
limYmin = -0.1*max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx));     

limYmaxMap = max( max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx)) ;
limYminMap = min( min(ATFObj(:,round(size(ATFObj,1)/2))),min(signalLambda0twoWingsfx)) ;

figure('Units', 'normalized', 'Position', [0.1 0.1 0.4 0.6])
    hAxis(1) = subplot(2,2,1)
        pcolor(fx*1E-6,fy*1E-6,  Illum_mat/max(max(Illum_mat)))
        rectangle('Position', [limfx-ATFbuffer-ATFbarlength, -limfx+ATFbuffer, ATFbarlength, ATFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
        xlim([-limfx limfx])
        ylim([-limfx limfx])
        shading flat    
        axis square
        xticks([])
        yticks([])
        set(gca, 'Fontsize', sizeFont)
        %c1 = colorbar
        %c1.Label.String = 'Amplitude'; 
        colormap(map)
        caxis([limYminMap limYmaxMap])
        set(gcf, 'Color', [1 1 1])
    
    hAxis(2) = subplot(2,2,2)
        pcolor(fx*1E-6,fy*1E-6, ATFObj)
        xlim([-limfx limfx])
        ylim([-limfx limfx])
        shading flat
        rectangle('Position', [limfx-ATFbuffer-ATFbarlength, -limfx+ATFbuffer, ATFbarlength, ATFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
        axis square
        xticks([])
        yticks([])
        set(gca, 'Fontsize', sizeFont)
        c1 = colorbar
        c1.Label.String = 'magnitude';   
        colormap(map)
        caxis([limYminMap limYmaxMap])
        set(gcf, 'Color', [1 1 1])
    
    hAxis(3) = subplot(2,2,[3 4])
        hold on
        area(fx*1E-6, signalLambda0twoWingsfx/max(signalLambda0twoWingsfx), 'Linewidth', 1.5, 'Facecolor', [1 0.8 0.8], 'Edgecolor', 'none');
        area(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Facecolor', [0.8 0.8 0.8], 'Edgecolor', 'none', 'FaceAlpha', 0.5);
        plot(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Color', [0 0 0]);
        plot(fx*1E-6, signalLambda0twoWingsfx/max(signalLambda0twoWingsfx), 'Linewidth', 1.5, 'Color', [1 0 0]);
        line([-limfx*1E6 limfx*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        
        xlim([-limfx limfx])
        ylim([-0.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx))) 1.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx)))])
        set(gca, 'Fontsize', sizeFont)
        xlabel('spatial frequency (1/\mum)');
        ylabel('magnitude')
        
        box on;
            
        pad = 0.1;
        padHor = 0.06;
        bottom = 0.4;
        
        sizeTotalPad = (2+bottom+2*pad);

        pos = get( hAxis(3), 'Position' )
        pos(2) = bottom/ sizeTotalPad ;                         % Shift down.
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(3), 'Position', pos ) 

        pos = get( hAxis(1), 'Position' )
        pos(1) = padHor;
        pos(2) = (1 + pad + bottom)/ sizeTotalPad ; 
        pos(3) = 0.45; % Shift down.
        pos(4) = 0.85/sizeTotalPad ;                        % Increase height.
        set( hAxis(1), 'Position', pos )

        pos = get( hAxis(2), 'Position' )
        pos(1) = 0.5-2*padHor;
        pos(2) = (1 + pad + bottom)/ sizeTotalPad ; 
        pos(3) = 0.45; % Shift down.
        pos(4) = 0.85/sizeTotalPad ;                        % Increase height.
        set( hAxis(2), 'Position', pos )  

        set(gcf, 'Color', [1 1 1])
    

        
%% Suppl Info plot 

limYmax = 1.1*max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx));
limYmin = -0.1*max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx));     

limYmaxMap = max( max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx)) ;
limYminMap = min( min(ATFObj(:,round(size(ATFObj,1)/2))),min(signalLambda0twoWingsfx)) ;

sizeFont = 18;
figure('Units', 'normalized', 'Position', [0.1 0.1 0.5 0.6])
      set(gcf, 'Color', [1 1 1])
      
%     %rectangle('Position',[0 0 1 1],'FaceColor',[1 1 0]);
      hAxis(5) = axes('Position',[0 0 1 0.5])
      set(gca,'Color','w')
      xticks([])
      yticks([])
      box off
      rectangle('Position',[0 0 1 1],'FaceColor',[1 0.9 0.9], 'EdgeColor', 'none');
      set(gca,'Visible','off')
      

    hAxis(1) = axes;%subplot(2,2,1)
        pcolor(fx*1E-6,fy*1E-6, ATFObj)
        xlim([-limfx limfx])
        ylim([-limfx limfx])
        shading flat
        %rectangle('Position', [limfx-ATFbuffer-ATFbarlength, -limfx+ATFbuffer, ATFbarlength, ATFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
        xlabel('u (1/µm)')
        ylabel('v (1/µm)')
        axis square
        %xticks([])
        %yticks([])
        set(gca, 'Fontsize', sizeFont)
        c1 = colorbar
        c1.Label.String = 'magnitude';   
        colormap(hAxis(1), gray)
        c1.FontSize = sizeFont;
        caxis([limYminMap limYmaxMap])
       
    
    
    
    hAxis(2) = axes;%subplot(2,2,2)
        hold on
        plot(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Color', [0 0 0]);
        line([-limfx*1E6 limfx*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        
        xlim([-limfx limfx])
        ylim([-0.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx))) 1.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx)))])
        set(gca, 'Fontsize', sizeFont)
        xlabel('cycles 1/\mum');
        ylabel('magnitude')    
        box on;
        axis square;
   
        
    hAxis(3) = axes;%subplot(2,2,3)
        pcolor(fx*1E-6,fy*1E-6,  Illum_mat/max(max(Illum_mat)))
        %rectangle('Position', [limfx-ATFbuffer-ATFbarlength, -limfx+ATFbuffer, ATFbarlength, ATFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
        xlim([-limfx limfx])
        ylim([-limfx limfx])
        shading flat    
        axis square
        %xticks([])
        %yticks([])
        
        c1 = colorbar
        c1.Label.String = 'magnitude'; 
        c1.FontSize = sizeFont;
        colormap(hAxis(3), map)
        set(gca, 'Fontsize', sizeFont)
        caxis([limYminMap limYmaxMap])
        set(gcf, 'Color', [1 1 1])
        xlabel('u (1/µm)')
        ylabel('v (1/µm)')
    
    hAxis(4) = axes;%subplot(2,2,4)
        hold on
        %area(fx*1E-6, signalLambda0twoWingsfx/max(signalLambda0twoWingsfx), 'Linewidth', 1.5, 'Facecolor', [1 0.8 0.8], 'Edgecolor', 'none');
        %area(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Facecolor', [0.8 0.8 0.8], 'Edgecolor', 'none', 'FaceAlpha', 0.5);
        %plot(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Color', [0 0 0]);
        plot(fx*1E-6, signalLambda0twoWingsfx/max(signalLambda0twoWingsfx), 'Linewidth', 1.5, 'Color', [1 0 0]);
        line([-limfx*1E6 limfx*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        
        xlim([-limfx limfx])
        ylim([-0.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx))) 1.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx)))])
        set(gca, 'Fontsize', sizeFont)
        xlabel('u (1/\mum)');
        ylabel('magnitude')
        box on;
        axis square;
            
        pad = 0.1;
        padHor = 0.06;
        bottom = 0.4;
        
        sizeTotalPad = (2+bottom+2*pad);

        pos = get( hAxis(1), 'Position' )
        pos(1) = 0.05;
        %pos(2) = (1 + pad + bottom) / sizeTotalPad ; 
        pos(3) = 0.4; % Shift down.
        pos(4) = 0.85/sizeTotalPad ;  
        pos(2) = (2-pos(4)*sizeTotalPad + pad + bottom) / sizeTotalPad ;
        set( hAxis(1), 'Position', pos )
        
        pos = get( hAxis(2), 'Position' )
        pos(1) = 0.65-2*padHor;
        pos(3) = 0.4;
        pos(4) = 0.85/sizeTotalPad ;  
        pos(2) = (2-pos(4)*sizeTotalPad + pad + bottom) / sizeTotalPad ;                  
        set( hAxis(2), 'Position', pos )
        
        pos = get( hAxis(3), 'Position' )
        pos(1) = 0.05;
        pos(2) = 0.85*bottom/ sizeTotalPad ;  
        pos(3) = 0.4; % Shift down.
        pos(4) = 0.85/sizeTotalPad ;  
        set( hAxis(3), 'Position', pos )
        
        
        pos = get( hAxis(4), 'Position' )
        pos(1) = 0.65-2*padHor;
        pos(2) = 0.85*bottom/ sizeTotalPad ;  
        pos(3) = 0.4; % Shift down.
        pos(4) = 0.85/sizeTotalPad ;  
        set( hAxis(4), 'Position', pos )
        
%         pos = get( hAxis(3), 'Position' )
%         pos(2) = bottom/ sizeTotalPad ;                         % Shift down.
%         pos(4) = 1/sizeTotalPad ;                        % Increase height.
%         set( hAxis(3), 'Position', pos ) 
% 
%         

%         pos = get( hAxis(4), 'Position' )
%         pos(1) = 0.5-2*padHor;
%         pos(2) = (1 + pad + bottom)/ sizeTotalPad ; 
%         pos(3) = 0.45; % Shift down.
%         pos(4) = 0.85/sizeTotalPad ;                        % Increase height.
%         set( hAxis(4), 'Position', pos )  

       
    

 %%
 
limYmax = 1.1 * max( max(max(ASFrealLightCorNorm)), max(max(ASFObjTheoNorm)) );
limYmin = min(min(min(ASFrealLightCorNorm)), min(min(ASFObjTheoNorm))) - 0.1/1.1*limYmax;     

limYmaxMap = max( max(max(ASFrealLightCorNorm)), max(max(ASFObjTheoNorm)) );
limYminMap = min( min(min(ASFrealLightCorNorm)), min(min(ASFObjTheoNorm)) );

 
 figure('Units', 'normalized', 'Position', [0.1 0.1 0.5 0.6])
      set(gcf, 'Color', [1 1 1])
      
%     %rectangle('Position',[0 0 1 1],'FaceColor',[1 1 0]);
      hAxis(5) = axes('Position',[0 0 1 0.5])
      set(gca,'Color','w')
      xticks([])
      yticks([])
      box off
      rectangle('Position',[0 0 1 1],'FaceColor',[1 0.9 0.9], 'EdgeColor', 'none');
      set(gca,'Visible','off')
      
      hAxis(1) = axes;%subplot(2,2,1)
        pcolor(xASFObj*1E6, yASFObj*1E6, ASFObjTheoNorm)
        
        xlim([-axlim*1E6 axlim*1E6])
        ylim([-axlim*1E6 axlim*1E6])
        shading flat    
        xlabel('x_2 - x_1 (\mum)');
        ylabel('y_2 - y_1 (\mum)');
        axis square

        set(gca, 'Fontsize', sizeFont)
        c1 = colorbar
        c1.Label.String = ('j_{bf}(x_2 - x_1,y_2 - y_1)');   
        colormap(hAxis(1), gray)
        c1.FontSize = sizeFont;
        caxis([limYminMap limYmaxMap])
        %set(gcf, 'Color', [1 1 1])
    
    hAxis(2) = axes; %subplot(2,2,2)
        hold on
        plot(xASFObj*1E6, ASFObjTheoNorm(:,round(size(ASFObjTheo,2)/2)), 'Linewidth', 1.5, 'Color', [0 0 0]) 
        line([-axlim*1E6 axlim*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        xlim([-axlim*1E6 axlim*1E6])
        ylim([limYmin limYmax])
        
        %xticks([-2 -1.5 -1, -0.5, 0, 0.5, 1 1.5 2])
        set(gca, 'Fontsize', sizeFont)
        xlabel('x_2 - x_1 (\mum)');
        ylabel('j_{bf}(x_2 - x_1,0)','visible','on');
        
        box on; 
        axis square;    
    
        
    hAxis(3) = axes;%subplot(2,2,3)
        pcolor(xASFObj*1E6, yASFObj*1E6, ASFrealLightCorNorm)
        xlim([-axlim*1E6 axlim*1E6])
        ylim([-axlim*1E6 axlim*1E6])
        shading flat    
        axis square
        
        set(gca, 'Fontsize', sizeFont)
        xlabel('x_2 - x_1 (\mum)');
        ylabel('y_2 - y_1 (\mum)');
        colormap(hAxis(3), map)
        caxis([limYminMap limYmaxMap])
        set(gca, 'Fontsize', sizeFont)
        c1 = colorbar
        c1.Label.String = ('j_{SEDFI}(x_2 - x_1,y_2 - y_1)')
        caxis([limYminMap limYmaxMap])    
    
    
    
        
    hAxis(4) = axes;%subplot(2,2,4)
        hold on
        %area(fx*1E-6, signalLambda0twoWingsfx/max(signalLambda0twoWingsfx), 'Linewidth', 1.5, 'Facecolor', [1 0.8 0.8], 'Edgecolor', 'none');
        %area(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Facecolor', [0.8 0.8 0.8], 'Edgecolor', 'none', 'FaceAlpha', 0.5);
        %plot(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Color', [0 0 0]);
        plot(xASFObj*1E6, ASFrealLightCorNorm(:,round(size(ASFrealLightCor,2)/2)), 'Linewidth', 1.5, 'Color', [1 0 0])
        line([-axlim*1E6 axlim*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        xlim([-axlim*1E6 axlim*1E6])
        ylim([limYmin limYmax])
        %xticks([-2 -1.5 -1, -0.5, 0, 0.5, 1 1.5 2])
        set(gca, 'Fontsize', sizeFont)
        xlabel('x_2 - x_1 (\mum)');
        ylabel('j_{SEDFI}(x_2 - x_1,0)','visible','on');
        
        box on; 
        axis square;
         
        
        pad = 0.1;
        padHor = 0.06;
        bottom = 0.4;
        
        sizeTotalPad = (2+bottom+2*pad);

        pos = get( hAxis(1), 'Position' )
        pos(1) = 0.05;
        %pos(2) = (1 + pad + bottom) / sizeTotalPad ; 
        pos(3) = 0.4; % Shift down.
        pos(4) = 0.85/sizeTotalPad ;  
        pos(2) = (2-pos(4)*sizeTotalPad + pad + bottom) / sizeTotalPad ;
        set( hAxis(1), 'Position', pos )
        
        pos = get( hAxis(2), 'Position' )
        pos(1) = 0.73-2*padHor;
        pos(3) = 0.4;
        pos(4) = 0.85/sizeTotalPad ;  
        pos(2) = (2-pos(4)*sizeTotalPad + pad + bottom) / sizeTotalPad ;                  
        set( hAxis(2), 'Position', pos )
        
        pos = get( hAxis(3), 'Position' )
        pos(1) = 0.05;
        pos(2) = 0.9*bottom/ sizeTotalPad ;  
        pos(3) = 0.4; % Shift down.
        pos(4) = 0.85/sizeTotalPad ;  
        set( hAxis(3), 'Position', pos )
        
        
        pos = get( hAxis(4), 'Position' )
        pos(1) = 0.73-2*padHor;
        pos(2) = 0.9*bottom/ sizeTotalPad ;  
        pos(3) = 0.4; % Shift down.
        pos(4) = 0.85/sizeTotalPad ;  
        set( hAxis(4), 'Position', pos )
        
%         pos = get( hAxis(3), 'Position' )
%         pos(2) = bottom/ sizeTotalPad ;                         % Shift down.
%         pos(4) = 1/sizeTotalPad ;                        % Increase height.
%         set( hAxis(3), 'Position', pos ) 
% 
%         

%         pos = get( hAxis(4), 'Position' )
%         pos(1) = 0.5-2*padHor;
%         pos(2) = (1 + pad + bottom)/ sizeTotalPad ; 
%         pos(3) = 0.45; % Shift down.
%         pos(4) = 0.85/sizeTotalPad ;                        % Increase height.
%         set( hAxis(4), 'Position', pos )  

    
    

       
        
        
        
   %% 
plotlimits = zeros(4,2); 
%%

%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_ls_real_res50nm.mat', 'x', 'y', 'IntIm');  
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_ls_real_res50nm_bf_diameter1000nm.mat', 'x', 'y', 'IntIm');  
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.4_ls_real_res50nm_df_diameter1000nm.mat', 'x', 'y', 'IntIm');  
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.37_ls_real_res50nm_df_diameter1000nm.mat'); 
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.43_ls_real_res50nm_df_diameter1000nm.mat'); 
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.46_ls_real_res50nm_df_diameter1000nm.mat'); 
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.58_ls_real_res50nm_df_diameter1000nm.mat'); 
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.85_ls_real_res50nm_df_diameter1000nm.mat'); 
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.7_ls_real_res50nm_df_diameter1000nm.mat'); 

%DF = 1; BF = 0;

%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.37_ls_real_res50nm_bf_diameter1000nm.mat', 'x', 'y', 'IntIm');  
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.4_ls_real_res50nm_bf_diameter1000nm.mat', 'x', 'y', 'IntIm');  
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.43_ls_real_res50nm_bf_diameter1000nm.mat', 'x', 'y', 'IntIm');  
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.46_ls_real_res50nm_bf_diameter1000nm.mat', 'x', 'y', 'IntIm');  
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.58_ls_real_res50nm_bf_diameter1000nm.mat', 'x', 'y', 'IntIm');  
%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.85_ls_real_res50nm_bf_diameter1000nm.mat', 'x', 'y', 'IntIm');  
load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.7_ls_real_res50nm_bf_diameter1000nm.mat', 'x', 'y', 'IntIm');  

DF = 0; BF = 1;

%load('/Users/mathias/Documents/MATLAB/Ceciles_paper/Jul2019/Oven/n120_nobj1.37_ls_real_res25nm_df_diameter1000nn.mat', 'x', 'y', 'IntIm');  

axlim = 1.3E-6;
fac = 2;
figure()
    pcolor(x(min(find(abs(x+axlim)<fac*min(abs(x)))):max(find(abs(x-axlim)<fac*min(abs(x))))), y(min(find(abs(y+axlim)<fac*min(abs(y)))):max(find(abs(y-axlim)<fac*min(abs(y))))), abs(IntIm(min(find(abs(x+axlim)<fac*min(abs(x)))):max(find(abs(x-axlim)<fac*min(abs(x)))), min(find(abs(y+axlim)<fac*min(abs(y)))):max(find(abs(y-axlim)<fac*min(abs(y))))))/max(max(abs(IntIm))))
    rectangle('Position', [axlim-ASFbuffer-ASFbarlength, -axlim+ASFbuffer, ASFbarlength, ASFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
   
    xlim([-axlim axlim])
    ylim([-axlim axlim])
    shading interp    
    axis square
    xticks([])
    yticks([])
    set(gca, 'Fontsize', sizeFont)
    c1 = colorbar
    c1.Label.String = 'Amplitude'; 
    colormap(map)
    set(gcf, 'Color', [1 1 1])

figure()
    plot(x(min(find(abs(x+axlim)<min(abs(x)))):max(find(abs(x-axlim)<min(abs(x))))), ...
           abs(IntIm(min(find(abs(x+axlim)<min(abs(x)))):max(find(abs(x-axlim)<min(abs(x)))), floor(size(IntIm,2)/2))/max(max(abs(IntIm)))));

       
mfactor = 1.005       

plotlimits(2,1+(DF~=1)) = max(max(abs(IntIm(min(find(abs(x+axlim)<fac*min(abs(x)))):max(find(abs(x-axlim)<fac*min(abs(x)))), floor(size(IntIm,2)/2))/max(max(abs(IntIm))))));
plotlimits(1,1+(DF~=1)) = min(min(abs(IntIm(min(find(abs(x+axlim)<fac*min(abs(x)))):max(find(abs(x-axlim)<fac*min(abs(x)))), floor(size(IntIm,2)/2))/max(max(abs(IntIm))))));

plotlimits(4,1+(DF~=1)) = mfactor*max(abs(IntIm(:,round(size(IntIm,1)/2)))/max(max(abs(IntIm))));
plotlimits(3,1+(DF~=1)) = min(min(abs(IntIm(min(find(abs(x+axlim)<fac*min(abs(x)))):max(find(abs(x-axlim)<fac*min(abs(x)))), floor(size(IntIm,2)/2))/max(max(abs(IntIm))))))-(mfactor-1)/mfactor*limYmax;

ASFplotCrop*max(xASFObj);

fraction = 0.01;
[Contrast, dContrast, nNumforAv, minMaxGrDF, maxMinGrDF] = findContrast(IntIm(min(find(abs(x+axlim)<fac*min(abs(x)))):max(find(abs(x-axlim)<fac*min(abs(x)))), round(size(IntIm,2)/2))/max(max(abs(IntIm))),fraction)



%%


figure('Units', 'normalized', 'Position', [0.1 0.1 0.3 0.6])
    hAxis(1) = subplot(2,1,1)
        pcolor(x, y, abs(IntIm)/max(max(abs(IntIm))))
        rectangle('Position', [axlim-ASFbuffer-ASFbarlength, -axlim+ASFbuffer, ASFbarlength, ASFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
        xlim([-axlim axlim])
        ylim([-axlim axlim])
        shading flat    
        axis square
        xticks([])
        yticks([])
        set(gca, 'Fontsize', sizeFont)
        c1 = colorbar
        c1.Label.String = 'amplitude';   
        if(DF==1) 
            colormap(map)
        else 
            colormap gray
        end
        caxis([min(plotlimits(1,:)) max(plotlimits(2,:))])
        set(gcf, 'Color', [1 1 1])
    
    hAxis(2) = subplot(2,1,2)
        hold on
        %plot(x*1E6, IntIm(:,round(size(IntIm,2)/2))/max(max(abs(IntIm))), 'Linewidth', 1.5, 'Color', [1 0 0]) 
        plot(x*1E6, IntIm(:,round(size(IntIm,2)/2))/max(max(abs(IntIm))), 'Linewidth', 1.5, 'Color', [(DF==1) 0 0]) 
        line([min(x*1E6) max(x*1E6)], [maxMinGrDF maxMinGrDF], 'Linestyle', ':', 'Linewidth', 1.5, 'Color', [0 0 0]);
        line([min(x*1E6) max(x*1E6)], [minMaxGrDF minMaxGrDF], 'Linestyle', ':', 'Linewidth', 1.5, 'Color', [0 0 0]);
        
        
        line([-axlim*1E6 axlim*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        xlim([-axlim*1E6 axlim*1E6])
        ylim([min(plotlimits(3,:)) max(plotlimits(4,:))])
        set(gca, 'Fontsize', sizeFont)
        xlabel('distance (\mum)');
        ylabel('amplitude','visible','on'); 
        xticks([-1, -0.5, 0, 0.5, 1])
        
        box on;
            
        pad = 0.1;
        padHor = 0.257;
        bottom = 0.4;
        
        sizeTotalPad = (2+bottom+2*pad);

        pos = get( hAxis(2), 'Position' )
        pos(1) = padHor ;         
        pos(2) = bottom/ sizeTotalPad ;
        pos(3) = pos(3) - padHor -0.038;% Shift down.
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(2), 'Position', pos ) 

        pos = get( hAxis(1), 'Position' )
        pos(1) = 0;
        pos(2) = (1 + pad + bottom)/ sizeTotalPad ; 
        pos(3) = 1; 
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(1), 'Position', pos )

        set(gcf, 'Color', [1 1 1])       

        
%%

nExp = size(x,2);
dxExp = x(2) - x(1);
fxExp = 1/dxExp*((-1/2+1/(2*nExp)):(1/nExp):(1/2-1/(2*nExp)));

IntImFFT = fftshift(fft2(IntIm));
IntImFFT(61,61) = 0;
figure()
    pcolor(fxExp, fxExp, (abs(IntImFFT)))
    shading flat
    colormap (map)
%%        
IntImSlice =  IntIm(:,round(size(IntIm,2)/2))/max(max(abs(IntIm)));
IntImSliceFFT = fftshift(fft(IntImSlice));
figure()
    plot(fxExp, real(IntImSliceFFT))
   
%% Colloid figures 
  sizeFont = 24;
DistperPix = 0.025E-6; %0.025µm
n1 = 1.37; % refractive index of object
lambda0 = 630E-9; % 630nm red

refrInd = 1.33; %refr. ind. of medium

% Alternative image to start with
n = 101; %120 takes 27.5h
x = DistperPix*(-n/2+1/2:1:n/2-1/2);
y = x';

radius = 20; %should correspond to a diameter of 1µm  

ObjectRaw = 2*sqrt((sqrt(x.^2 + y.^2) < radius*DistperPix).*((radius*DistperPix)^2 - (x.^2 + y.^2)));
Object = exp(1i*2*pi/lambda0 * (n1-refrInd).*ObjectRaw);


figure()
  pcolor(x, y, unwrap(angle(Object)))
  rectangle('Position', [axlim-ASFbuffer-ASFbarlength, -axlim+ASFbuffer, ASFbarlength, ASFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
  xlim([min(x) max(x)])
  ylim([min(y) max(y)])
        shading flat    
        axis square
        xticks([])
        yticks([])
        set(gca, 'Fontsize', sizeFont)
        %c1 = colorbar
        %c1.Label.String = 'Amplitude'; 
        colormap(map)
        set(gcf, 'Color', [1 1 1])


colormap gray
cmap = colormap;    
cmaptop = cmap;
cmaptop(:,1) = 1;
cmaptop(:,2) = 1;

cmapbottom = cmap;
cmapbottom(:,3) = 0;
map = [cmapbottom; cmaptop];
map = 0.95*map;


limYmax = 1.2*max(max(unwrap(angle((Object)))))/2/pi;
limYmin = min(min(unwrap(angle((Object))))) - 0.1/1.1*limYmax;     

limYmaxMap = max(max(unwrap(angle(Object))))/2/pi;
limYminMap = min(min(unwrap(angle(Object))))/2/pi;

axlim = max(x);

figure('Units', 'normalized', 'Position', [0.1 0.1 0.26 0.6])
    hAxis(1) = subplot(2,1,1)
        pcolor(x, y, unwrap(angle(Object))/2/pi)
        rectangle('Position', [axlim-ASFbuffer-ASFbarlength, -axlim+ASFbuffer, ASFbarlength, ASFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
        xlim([-axlim axlim])
        ylim([-axlim axlim])
        shading flat    
        axis square
        xticks([])
        yticks([])
        set(gca, 'Fontsize', sizeFont)
        colormap(map)
        set(gcf, 'Color', [1 1 1])
        c1 = colorbar
        c1.Label.String = 'phase / 2\pi';   
        colormap(map)
        caxis([limYminMap limYmaxMap])
        set(gcf, 'Color', [1 1 1])
    
    hAxis(2) = subplot(2,1,2)
        hold on
        area(x*1E6, abs(Object(:, floor(size(Object,1)/2))), 'Linewidth', 1.5, 'Facecolor', [0.8 0.8 0.8], 'Edgecolor', 'none', 'FaceAlpha', 0.5) 
        plot(x*1E6, abs(Object(:, floor(size(Object,1)/2))), 'Linewidth', 1.5, 'Color', [0 0 0]) 
        
        xlim([-axlim*1E6 axlim*1E6])
        ylim([-0.1 1.1])
        yticks([0 0.25 0.5 0.75 1])
        ylabel('magnitude','visible','on');
        yyaxis right
        area(x*1E6, unwrap(angle(Object(:, floor(size(Object,1)/2))))/2/pi, 'Facecolor', [0.8 0.8 0], 'Edgecolor', 'none', 'FaceAlpha', 0.2);      
        plot(x*1E6, unwrap(angle(Object(:, floor(size(Object,1)/2))))/2/pi, 'Linewidth', 1.5, 'Color', [0.6 0.6 0])
        line([-axlim*1E6 axlim*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        ylim([limYmin limYmax])
        
        set(gca, 'Fontsize', sizeFont)
        xlabel('distance (\mum)');
        ylabel('phase / 2\pi','visible','on');
        xticks([-1, -0.5, 0, 0.5, 1])
        
        set(hAxis(2), 'XColor', 'k', 'YColor', [0.6 0.6 0]);
        
        box on;
            
        pad = 0.1;
        padHor = 0.257;
        bottom = 0.4;
        
        sizeTotalPad = (2+bottom+2*pad);

        pos = get( hAxis(2), 'Position' )
        pos(1) = padHor;
        pos(2) = bottom/ sizeTotalPad ;                         % Shift down.
        pos(3) = pos(3) - padHor -0.038;% Shift down.
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(2), 'Position', pos ) 

        pos = get( hAxis(1), 'Position' )
        pos(1) = 0;
        pos(2) = (1 + pad + bottom)/ sizeTotalPad ; 
        pos(3) = 1; 
        pos(4) = 0.85/sizeTotalPad ;                        % Increase height.
        set( hAxis(1), 'Position', pos )

        set(gcf, 'Color', [1 1 1])
        
        
        
%% Microscope ASF ATF 

colormap gray
cmap = colormap;    
cmaptop = cmap;
cmaptop(:,3) = 1;

cmapbottom = cmap;
cmapbottom(:,1) = 0;
cmapbottom(:,2) = 0;
map = [cmapbottom; cmaptop];


% ASF 
limYmax = 1.1*max(ASFObjTheoNorm(:,round(size(ASFObjTheo,1)/2)));
limYmin = min(ASFObjTheoNorm(:,round(size(ASFObjTheo,1)/2)))-0.1/1.1*limYmax;

limYmaxMap = max(max(ASFObjTheoNorm));
limYminMap = min(min(ASFObjTheoNorm));

ASFplotCrop*max(xASFObj);

figure('Units', 'normalized', 'Position', [0.1 0.1 0.26 0.6])
    hAxis(1) = subplot(2,1,1)
        pcolor(xASFObj, yASFObj, ASFObjTheoNorm)
        rectangle('Position', [axlim-ASFbuffer-ASFbarlength, -axlim+ASFbuffer, ASFbarlength, ASFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
        xlim([-axlim axlim])
        ylim([-axlim axlim])
        shading flat    
        axis square
        xticks([])
        yticks([])
        set(gca, 'Fontsize', sizeFont)
        c1 = colorbar
        c1.Label.String = 'amplitude';   
        colormap (map)
        caxis([limYminMap limYmaxMap])
        set(gcf, 'Color', [1 1 1])
    
    hAxis(2) = subplot(2,1,2)
        hold on
        area(xASFObj*1E6, ASFObjTheoNorm(:,round(size(ASFObjTheo,2)/2)), 'Linewidth', 1.5, 'Facecolor', [0.8 0.8 1], 'Edgecolor', 'none')
        plot(xASFObj*1E6, ASFObjTheoNorm(:,round(size(ASFObjTheo,2)/2)), 'Linewidth', 1.5, 'Color', [0 0 1]) 
                
        line([-axlim*1E6 axlim*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        xlim([-axlim*1E6 axlim*1E6])
        ylim([limYmin limYmax])
        set(gca, 'Fontsize', sizeFont)
        xlabel('distance (\mum)');
        ylabel('amplitude','visible','on');
        xticks([-1, -0.5, 0, 0.5, 1])
        
        box on;
            
        pad = 0.1;
        padHor = 0.257;
        bottom = 0.4;
        
        sizeTotalPad = (2+bottom+2*pad);

        pos = get( hAxis(2), 'Position' )
        pos(1) = padHor;
        pos(2) = bottom/ sizeTotalPad ;                         % Shift down.
        pos(3) = pos(3) - padHor -0.038;% Shift down.
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(2), 'Position', pos ) 

        pos = get( hAxis(1), 'Position' )
        pos(1) = 0;
        pos(2) = (1 + pad + bottom)/ sizeTotalPad ; 
        pos(3) = 1; 
        pos(4) = 0.85/sizeTotalPad ;                        % Increase height.
        set( hAxis(1), 'Position', pos )

        set(gcf, 'Color', [1 1 1])
  
        
% ATF
limYmax = 1.1*max(ATFObj(:,round(size(ATFObj,1)/2)));
limYmin = min(ATFObj(:,round(size(ATFObj,1)/2)))-0.1/1.1*limYmax;

limYmaxMap = max(ATFObj(:,round(size(ATFObj,1)/2)));
limYminMap = min(ATFObj(:,round(size(ATFObj,1)/2))); 


figure('Units', 'normalized', 'Position', [0.1 0.1 0.25 0.6])
    hAxis(1) = subplot(2,1,1)
        pcolor(fx*1E-6,fy*1E-6, ATFObj)
        xlim([-limfx limfx])
        ylim([-limfx limfx])
        shading flat
        rectangle('Position', [limfx-ATFbuffer-ATFbarlength, -limfx+ATFbuffer, ATFbarlength, ATFbheight], 'Edgecolor', 'none', 'Facecolor', [1 1 1])
        axis square
        xticks([])
        yticks([])
        set(gca, 'Fontsize', sizeFont)
        c1 = colorbar
        c1.Label.String = 'magnitude';   
        colormap gray
        caxis([limYminMap limYmaxMap])
        set(gcf, 'Color', [1 1 1])
    
    hAxis(2) = subplot(2,1,2)
        hold on
        plot(fx*1E-6, ATFObj(:,round(size(ATFObj,1)/2)),'Linewidth', 1.5, 'Color', [0 0 0]);
        line([-limfx*1E6 limfx*1E6], [0 0], 'Linestyle', ':', 'Linewidth', 1, 'Color', [0 0 0])
        
        xlim([-limfx limfx])
        ylim([-0.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx))) 1.1*(max(max(ATFObj(:,round(size(ATFObj,1)/2))),max(signalLambda0twoWingsfx)))])
        set(gca, 'Fontsize', sizeFont)
        xlabel('cycles 1/\mum');
        ylabel('magnitude')
        
        box on;
            
        pad = 0.1;
        padHor = 0.12;
        bottom = 0.4;
        
        sizeTotalPad = (2+bottom+2*pad);

        pos = get( hAxis(2), 'Position' )
        pos(2) = bottom/ sizeTotalPad ;                         % Shift down.
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(2), 'Position', pos ) 

        pos = get( hAxis(1), 'Position' )
        pos(1) = 0;
        pos(2) = (1 + pad + bottom)/ sizeTotalPad ; 
        pos(3) = 1; 
        pos(4) = 1/sizeTotalPad ;                        % Increase height.
        set( hAxis(1), 'Position', pos )

        set(gcf, 'Color', [1 1 1])
          
        
    %%



% need to take care of the 1-pixel phase undulation in the Fourier transform of of the ATF in Illum_mat ! 
checkersY = -floor(size(ObjFFT,1)/2):1:floor(size(ObjFFT,1)/2);
checkersX = checkersY';
checkers = (2*mod(checkersX,2)-1)*(2*mod(checkersY,2)-1);
ObjFFTCor = checkers.*sign(imag(ObjFFT)).*angle(ObjFFT);
 

figure()
    pcolor(fSample, fSample, ObjFFTCor)
    shading flat    
    axis square
    xticks([])
    yticks([])
    set(gca, 'Fontsize', sizeFont)
    c1 = colorbar
    c1.Label.String = 'Amplitude'; 
    colormap gray
    set(gcf, 'Color', [1 1 1])

figure()
    pcolor(fSample, fSample, angle(ObjFFT))
    shading flat    
    axis square
    xticks([])
    yticks([])
    set(gca, 'Fontsize', sizeFont)
    c1 = colorbar
    c1.Label.String = 'Amplitude'; 
    colormap gray
    set(gcf, 'Color', [1 1 1])    