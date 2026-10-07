% MS and SLM PSF comparison and FWHM analysis.


clear
close all
clc

%% Parameters

checkFit = 0;
fitRange = 100;

magnification = 40*200/180;
pixelSize = 5.86/magnification;       % um/pixel

plotHalfRange = 20;                    % um
centerFitHalfRange = 7;                % um

%% Load MS and SLM measurements

D = load(fullfile(fileparts(mfilename('fullpath')), 'DarkFrame_100us.mat'));
S = load(fullfile(fileparts(mfilename('fullpath')), 'Focus_MS_20400um_120cores_100us.mat'));

img_MS = double(S.data_img) - double(D.DarkFrame);

D = load(fullfile(fileparts(mfilename('fullpath')), 'DarkFrame_500us.mat'));
S = load(fullfile(fileparts(mfilename('fullpath')), 'Focus_SLM_400um_120cores_500us.mat'));

img_SLM = double(S.data_img) - double(D.DarkFrame);

%% FWHM

[~,pos_MS,FWHMpix_MS] = ...
    GetStrehl_useMax_Generic(img_MS,checkFit,fitRange);

[~,pos_SLM,FWHMpix_SLM] = ...
    GetStrehl_useMax_Generic(img_SLM,checkFit,fitRange);

FWHM_MS = FWHMpix_MS * pixelSize;
FWHM_SLM = FWHMpix_SLM * pixelSize;

FWHM_MS_2p = FWHM_MS/sqrt(2);
FWHM_SLM_2p = FWHM_SLM/sqrt(2);

%% PSF images

img_MS = max(img_MS,0);
img_SLM = max(img_SLM,0);

img_MS_2p = img_MS.^2;
img_SLM_2p = img_SLM.^2;

clim_PSF = [0 max([img_MS(:); img_SLM(:)])];
clim_2p = [0 max([img_MS_2p(:); img_SLM_2p(:)])];

figure
imagesc(img_MS)
axis image off
colormap turbo
clim(clim_PSF)

figure
imagesc(img_MS_2p)
axis image off
colormap turbo
clim(clim_2p)

figure
imagesc(img_SLM)
axis image off
colormap turbo
clim(clim_PSF)

figure
imagesc(img_SLM_2p)
axis image off
colormap turbo
clim(clim_2p)

%% Horizontal PSF cuts

halfRange_px = round(plotHalfRange/pixelSize);

row = pos_MS(1);
col = pos_MS(2);
cols = max(1,col-halfRange_px):min(size(img_MS,2),col+halfRange_px);

x_MS = (cols-col)*pixelSize;
profile_MS = double(img_MS(row,cols));
profile_MS = (profile_MS-min(profile_MS)) / ...
    (max(profile_MS)-min(profile_MS));

[~,i0] = max(profile_MS);
fitMask = abs(x_MS-x_MS(i0)) <= centerFitHalfRange;
g = fit(x_MS(fitMask)',profile_MS(fitMask)','gauss1');
x_MS = x_MS-g.b1;

row = pos_SLM(1);
col = pos_SLM(2);
cols = max(1,col-halfRange_px):min(size(img_SLM,2),col+halfRange_px);

x_SLM = (cols-col)*pixelSize;
profile_SLM = double(img_SLM(row,cols));
profile_SLM = (profile_SLM-min(profile_SLM)) / ...
    (max(profile_SLM)-min(profile_SLM));

[~,i0] = max(profile_SLM);
fitMask = abs(x_SLM-x_SLM(i0)) <= centerFitHalfRange;
g = fit(x_SLM(fitMask)',profile_SLM(fitMask)','gauss1');
x_SLM = x_SLM-g.b1;

%% Excitation PSF resolution

figure
plot(x_MS,profile_MS,'x','MarkerSize',4,'LineWidth',0.8)
hold on
plot(x_SLM,profile_SLM,'-','LineWidth',1.4)
hold off

xlabel('X (\mum)')
ylabel('Normalized intensity')
title('Excitation PSF resolution')
legend(sprintf('MS, FWHM %.2f \\mum',FWHM_MS), ...
       sprintf('SLM, FWHM %.2f \\mum',FWHM_SLM), ...
       'Location','northeast','Box','off')
xlim([-20 20])
ylim([0 1])
grid on
box on

%% Two-photon signal resolution

profile_MS_2p = profile_MS.^2;
profile_SLM_2p = profile_SLM.^2;

profile_MS_2p = profile_MS_2p/max(profile_MS_2p);
profile_SLM_2p = profile_SLM_2p/max(profile_SLM_2p);

figure
plot(x_MS,profile_MS_2p,'x','MarkerSize',4,'LineWidth',0.8)
hold on
plot(x_SLM,profile_SLM_2p,'-','LineWidth',1.4)
hold off

xlabel('X (\mum)')
ylabel('Normalized I^2')
title('Two-photon signal resolution')
legend(sprintf('MS, FWHM %.2f \\mum',FWHM_MS_2p), ...
       sprintf('SLM, FWHM %.2f \\mum',FWHM_SLM_2p), ...
       'Location','northeast','Box','off')
xlim([-20 20])
ylim([0 1])
grid on
box on
