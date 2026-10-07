%% Post-processing of TPEF endoscopic images from PMT counts
% Input: PMTimg (2D matrix, counts per pixel)
% Output:
%   Fig 1: raw vs processed (not normalized)
%   Fig 2: raw, normalized only, no processing applied, transparent PNG export
%   Fig 3: raw normalized only, transparent PNG export
%   Fig 4: raw image with signal/background ROI contours (SNR/contrast analysis)
%   Fig 5: raw image with signal/background ROI contours, transparent PNG export
%
% Requires Image Processing Toolbox for some optional functions (drawpolygon,
% imgaussfilt, medfilt2, imopen)
% 
% % Fluo beads
% tmp = load('Beads_300us_5x5mRad.mat');
% fileName = 'FluorescentBeads';
% FOVx_um = 50;

% % Live Cell measurement
% tmp = load('LiveCell_10ms_6x6mRad.mat')
% fileName = 'LiveCell';
% FOVx_um = 60;

%Mouse Brain
tmp = load('BrainSlice_50ms_7x7mRad.mat')
fileName = 'MouseBrain';
FOVx_um = 70;

PMTimg = tmp.PMTimage;
PMTimg = PMTimg(2:end,2:end);
PMTimg = PMTimg ./ max(PMTimg,[],'all');
% PMTimg = PMTimg./max(PMTimg,[],'all');

clear tmp

% clearvars -except PMTimg
close all
clc

%% ----------------------- USER PARAMETERS -----------------------
SaveFigures = false;
% SaveFigures = false;

scaleBar_um = 10;       % e.g. 10 or 20

outFolder = ''; mkdir(outFolder);
outName = strcat(outFolder,fileName);

% ---------- Processing switches ----------
% doBaselineCorrection = false;
doBaselineCorrection = true;
doBorderRemoval      = false;
% doBorderRemoval      = true;
doFourierFilter      = true;
% doSpatialSmoothing   = false;
doSpatialSmoothing   = true;
doNormalizeDisplay   = true;

% ---------- SNR / contrast analysis (RAW image) ----------
doSNRContrastAnalysis = true;   % interactive: draws polygon ROIs, blocks execution

% ---------- Baseline correction ----------
baselineMethod = 'percentile';   % 'percentile' or 'morphopen'
baselinePercentile = 10;
morphRadius = 2;

% ---------- Border / outline removal ----------
borderWidth = 1;
% borderMode  = 'subtractMedian';  % 'subtractMedian' or 'replaceBorder'
borderMode  = 'replaceBorder';  % 'subtractMedian' or 'replaceBorder'

% ---------- Fourier low-pass filter -----------------------------------

fftFilterType = 'gaussian';      % 'gaussian' or 'butterworth'
cutoffFeatureSize_um = 1.0;      % smallest feature size (um) to preserve; re-tune empirically
butterOrder   = 2;

% ---------- Optional spatial smoothing ----------
spatialFilterType = 'gaussian';  % 'gaussian' or 'median'
gaussSigma        = 0.6;
medianKernel      = [3 3];

% ---------- Display colormap ----------
nCmap = 256;
greenMap = [zeros(nCmap,1), linspace(0,1,nCmap)', zeros(nCmap,1)];
yellowMap = [linspace(0,1,nCmap)', linspace(0,1,nCmap)', zeros(nCmap,1)];

% ---------- Display normalization ----------
clipDisplayToPercentiles = true;
displayLowPct  = 0.5;
displayHighPct = 99.9;

% ---------- Scalebar settings ----------
% addScaleBar = true;
addScaleBar = false;

if addScaleBar == true
    outName = strcat(outName,'_Scalebar');
else
    outName = strcat(outName,'_NoScalebar');
end



% Scalebar appearance
scaleBarCorner = 'southwest';   % 'southwest','southeast','northwest','northeast'
% scaleBarCorner = 'northeast';   % 'southwest','southeast','northwest','northeast'
scaleBarLineWidth = 4;
scaleBarTextColor = 'w';
scaleBarFontSize = 11;
scaleBarMarginPx = 3;           % distance from border in pixels
scaleBarTextGapPx = 1;          % vertical gap between line and text

% ---------- Export ----------
if SaveFigures == true
    saveComparisonPNG = false;
    saveComparisonFIG = false;
    
    saveProcessedTransparentPNG = true;
    saveProcessedFIG = false;
    saveProcessedTIFF = false;
    saveRawTransparentPNG = true;
    saveRawFIG = true;

    saveSNRContrastPNG = true;
    saveSNRContrastFIG = true;
    saveSNRContrastMAT = true;
else 
    saveComparisonPNG = false;
    saveComparisonFIG = false;
    
    saveProcessedTransparentPNG = false;
    saveProcessedFIG = false;
    saveProcessedTIFF = false;
    saveRawTransparentPNG = false;
    saveRawFIG = false;

    saveSNRContrastPNG = false;
    saveSNRContrastFIG = false;
    saveSNRContrastMAT = false;
    disp('NOT saving!')
end

%% ----------------------- INPUT CHECK -----------------------
if ~exist('PMTimg','var')
    error('Variable PMTimg is not found in the workspace.');
end

if ~ismatrix(PMTimg) || ~isnumeric(PMTimg)
    error('PMTimg must be a numeric 2D matrix.');
end

img0 = double(PMTimg);

if any(~isfinite(img0(:)))
    warning('PMTimg contains NaN or Inf values. Replacing with zero.');
    img0(~isfinite(img0)) = 0;
end

[nRows, nCols] = size(img0);

fprintf('Input image size: %d x %d pixels\n', nRows, nCols);

%% ----------------------- SCALEBAR CALCULATION -----------------------

if addScaleBar
    pixSizeX_um = FOVx_um / nCols;
    scaleBar_px = scaleBar_um / pixSizeX_um;

    % Round to nearest integer pixel
    scaleBar_px = round(scaleBar_px);

    if scaleBar_px < 2
        warning('Scalebar is shorter than 2 pixels. Increase FOV precision or use a larger scaleBar_um.');
    end

    if scaleBar_px > round(0.8*nCols)
        warning('Scalebar is too large relative to image width.');
    end
else
    scaleBar_px = [];   % unused downstream when addScaleBar is false
end

%% ----------------------- STEP 0: SNR / CONTRAST ANALYSIS (RAW IMAGE) -----

if doSNRContrastAnalysis
    snrContrastResults = analyzeSignalBackgroundROI(img0, greenMap, ...
        outName, saveSNRContrastPNG, saveSNRContrastFIG, saveSNRContrastMAT, ...
        addScaleBar, scaleBar_px, scaleBar_um, scaleBarCorner, scaleBarMarginPx, ...
        scaleBarLineWidth, scaleBarTextColor, scaleBarFontSize, scaleBarTextGapPx);

    fprintf('\n--- SNR / Contrast (raw image) ---\n');
    fprintf('  mean_signal      = %.2f\n', snrContrastResults.mean_signal);
    fprintf('  std_signal       = %.2f\n', snrContrastResults.std_signal);
    fprintf('  mean_background  = %.2f\n', snrContrastResults.mean_background);
    fprintf('  std_background   = %.2f\n', snrContrastResults.std_background);
    fprintf('  SNR (sig/std_sig)      = %.2f\n', snrContrastResults.SNR_signal_over_stdSignal);
    fprintf('  SNR (sig/std_bg)       = %.2f\n', snrContrastResults.SNR_signal_over_stdBackground);
    fprintf('  SNR corrected ((sig-bg)/std_bg), additional, not requested = %.2f\n', snrContrastResults.SNR_corrected);
    fprintf('  Contrast (sig/bg)      = %.2f\n\n', snrContrastResults.contrast);
end

%% ----------------------- STEP 1: BASELINE CORRECTION -----------------------
img1 = img0;

if doBaselineCorrection
    switch lower(baselineMethod)
        case 'percentile'
            bg = prctile(img1(:), baselinePercentile);
            img1 = img1 - bg;

        case 'morphopen'
            se = strel('disk', morphRadius, 0);
            bg = imopen(img1, se);
            img1 = img1 - bg;

        otherwise
            error('Unknown baselineMethod: %s', baselineMethod);
    end

    img1(img1 < 0) = 0;
end

%% ----------------------- STEP 2: BORDER / OUTLINE REMOVAL -----------------------
img2 = img1;

if doBorderRemoval
    bw = borderWidth;

    if 2*bw >= min(nRows, nCols)
        error('borderWidth is too large for this image size.');
    end

    borderMask = false(nRows, nCols);
    borderMask(1:bw,:) = true;
    borderMask(end-bw+1:end,:) = true;
    borderMask(:,1:bw) = true;
    borderMask(:,end-bw+1:end) = true;

    borderVals = img2(borderMask);
    borderMedian = median(borderVals(:));

    switch lower(borderMode)
        case 'subtractmedian'
            img2 = img2 - borderMedian;
            img2(img2 < 0) = 0;

        case 'replaceborder'
            interior = img2(bw+1:end-bw, bw+1:end-bw);
            fillVal = median(interior(:));
            img2(borderMask) = fillVal;

        otherwise
            error('Unknown borderMode: %s', borderMode);
    end
end

%% ----------------------- STEP 3: FOURIER FILTERING (fixed spatial frequency) --
img3 = img2;

if doFourierFilter
    F = fftshift(fft2(img3));

    pixSizeX_um = FOVx_um / nCols;
    pixSizeY_um = FOVx_um / nRows;

    % Frequency axes in cycles/um (physical units), not fractions of Nyquist
    fx_axis = (-floor(nCols/2):ceil(nCols/2)-1) / (nCols*pixSizeX_um);
    fy_axis = (-floor(nRows/2):ceil(nRows/2)-1) / (nRows*pixSizeY_um);
    [FX, FY] = meshgrid(fx_axis, fy_axis);
    FR = sqrt(FX.^2 + FY.^2);   % radial spatial frequency, cycles/um

    cutoffFreq = 1 / cutoffFeatureSize_um;   % cycles/um

    switch lower(fftFilterType)
        case 'gaussian'
            sigmaF = cutoffFreq;
            H = exp(-(FR.^2) / (2*sigmaF^2));

        case 'butterworth'
            D0 = cutoffFreq;
            H = 1 ./ (1 + (FR ./ D0).^(2*butterOrder));

        otherwise
            error('Unknown fftFilterType: %s', fftFilterType);
    end

    Ffilt = F .* H;
    % figure(99)
    % imagesc(fx_axis, fy_axis, real(Ffilt))
    % pbaspect([1 1 1])
    % title('Fourier domain')
    img3 = real(ifft2(ifftshift(Ffilt)));
    img3(img3 < 0) = 0;
end

%% ----------------------- STEP 4: OPTIONAL SPATIAL SMOOTHING -----------------------
img4 = img3;

if doSpatialSmoothing
    switch lower(spatialFilterType)
        case 'gaussian'
            img4 = imgaussfilt(img4, gaussSigma);

        case 'median'
            img4 = medfilt2(img4, medianKernel);

        otherwise
            error('Unknown spatialFilterType: %s', spatialFilterType);
    end
end

%% ----------------------- STEP 5: OUTPUT IMAGES -----------------------
imgProcessed = img4;   % absolute processed image, not normalized

imgDisp = imgProcessed;
if doNormalizeDisplay
    if clipDisplayToPercentiles
        lo = prctile(imgDisp(:), displayLowPct);
        hi = prctile(imgDisp(:), displayHighPct);
    else
        lo = min(imgDisp(:));
        hi = max(imgDisp(:));
    end

    if hi <= lo
        lo = min(imgDisp(:));
        hi = max(imgDisp(:));
    end

    if hi > lo
        imgDisp = (imgDisp - lo) / (hi - lo);
    else
        imgDisp = zeros(size(imgDisp));
    end

    imgDisp = max(0, min(1, imgDisp));
end

%% ----------------------- FIGURE 1: RAW VS PROCESSED -----------------------
fig1 = figure('Color','w','Name','Raw_vs_Processed');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

ax1 = nexttile;
imagesc(ax1, img0, [0 1])
axis image off
title('Raw')
colormap(ax1, greenMap)
% colormap(ax1, yellowMap)
colorbar

ax2 = nexttile;
imagesc(ax2, imgProcessed./max(imgProcessed,[],"all"), [0 1])
axis image off
title('Processed')
colormap(ax2, greenMap)
% colormap(ax2, yellowMap)
colorbar
hold(ax2,'on')

if addScaleBar
    drawScaleBar(ax2, nRows, nCols, scaleBar_px, scaleBar_um, ...
        scaleBarCorner, scaleBarMarginPx, scaleBarLineWidth, ...
        scaleBarTextColor, scaleBarFontSize, scaleBarTextGapPx);
end

% sgtitle(sprintf('TPEF image processing (%d x %d)', nRows, nCols))
drawnow

%% ----------------------- FIGURE 2: RAW, NORMALIZED ONLY (NO PROCESSING) -----------------------
% Transparent figure/background export. Plots img0 directly - NOT
% imgProcessed/imgDisp - so no baseline correction, Fourier filtering,
% spatial smoothing, or percentile display-stretching is applied here.
% img0 is already at max = 1 thanks to the top-of-script PMTimg
% normalization, so no further division is needed.
fig2 = figure('Color','none','Name','Raw_Normalized_Transparent');
ax3 = axes(fig2);
imagesc(ax3, img0, [0 1])
axis(ax3, 'image')
axis(ax3, 'off')
colormap(ax3, greenMap)
% colormap(ax3, yellowMap)
colorbar(ax3)
set(ax3, 'Color', 'none')
hold(ax3, 'on')

if addScaleBar
    drawScaleBar(ax3, nRows, nCols, scaleBar_px, scaleBar_um, ...
        scaleBarCorner, scaleBarMarginPx, scaleBarLineWidth, ...
        scaleBarTextColor, scaleBarFontSize, scaleBarTextGapPx);
end

drawnow

%% ----------------------- FIGURE 3: NORMALIZED RAW ONLY -----------------------
% Transparent figure/background export
fig3 = figure('Color','none','Name','Raw_Normalized_Transparent');
ax4 = axes(fig3);
imagesc(ax4, img0, [0 1])
axis(ax4, 'image')
axis(ax4, 'off')
colormap(ax4, greenMap)
% colormap(ax4, yellowMap)
set(ax4, 'Color', 'none')
hold(ax4, 'on')

if addScaleBar
    drawScaleBar(ax4, nRows, nCols, scaleBar_px, scaleBar_um, ...
        scaleBarCorner, scaleBarMarginPx, scaleBarLineWidth, ...
        scaleBarTextColor, scaleBarFontSize, scaleBarTextGapPx);
end

drawnow

%% ----------------------- SAVE OUTPUTS -----------------------
% Figure 1: comparison
if saveComparisonPNG
    exportgraphics(fig1, [outName '_comparison.png'], 'Resolution', 300);
end

if saveComparisonFIG
    savefig(fig1, [outName '_comparison.fig']);
end

% Figure 2: normalized processed only, transparent PNG
if saveProcessedTransparentPNG
    exportgraphics(ax3, [outName '_processed_transparent.png'], ...
        'BackgroundColor', 'none', 'Resolution', 300);
end

if saveProcessedFIG
    savefig(fig2, [outName '_processed.fig']);
end

% TIFF from normalized processed image
if saveProcessedTIFF
    img16 = uint16(65535 * imgDisp);
    imwrite(img16, [outName '_processed.tif'], 'tif', 'Compression', 'none');
end

% Figure 3: normalized raw only, transparent PNG
if saveRawTransparentPNG
    exportgraphics(ax4, [outName '_raw_transparent.png'], ...
        'BackgroundColor', 'none', 'Resolution', 300);
end

if saveRawFIG
    savefig(fig3, [outName '_raw.fig']);
end

fprintf('\nSaved files:\n');
if saveComparisonPNG
    fprintf('  %s\n', [outName '_comparison.png']);
end
if saveComparisonFIG
    fprintf('  %s\n', [outName '_comparison.fig']);
end
if saveProcessedTransparentPNG
    fprintf('  %s\n', [outName '_processed_transparent.png']);
end
if saveProcessedFIG
    fprintf('  %s\n', [outName '_processed.fig']);
end
if saveProcessedTIFF
    fprintf('  %s\n', [outName '_processed.tif']);
end
if saveRawTransparentPNG
    fprintf('  %s\n', [outName '_raw_transparent.png']);
end
if saveRawFIG
    fprintf('  %s\n', [outName '_raw.fig']);
end
if saveSNRContrastPNG
    fprintf('  %s\n', [outName '_ROI_overlay.png']);
    fprintf('  %s\n', [outName '_ROI_overlay_transparent.png']);
end
if saveSNRContrastFIG
    fprintf('  %s\n', [outName '_ROI_overlay.fig']);
end
if saveSNRContrastMAT
    fprintf('  %s\n', [outName '_SNR_contrast_results.mat']);
end

%% ----------------------- LOCAL FUNCTIONS -----------------------
function drawScaleBar(ax, nRows, nCols, scaleBar_px, scaleBar_um, ...
    corner, marginPx, lineWidth, textColor, fontSize, textGapPx)

    switch lower(corner)
        case 'southwest'
            x1 = 1 + marginPx;
            x2 = x1 + scaleBar_px;
            y  = nRows - marginPx;

            textX = (x1 + x2)/2;
            textY = y - textGapPx;

            vAlign = 'bottom';

        case 'southeast'
            x2 = nCols - marginPx;
            x1 = x2 - scaleBar_px;
            y  = nRows - marginPx;

            textX = (x1 + x2)/2;
            textY = y - textGapPx;

            vAlign = 'bottom';

        case 'northwest'
            x1 = 1 + marginPx;
            x2 = x1 + scaleBar_px;
            y  = 1 + marginPx;

            textX = (x1 + x2)/2;
            textY = y + textGapPx;

            vAlign = 'top';

        case 'northeast'
            x2 = nCols - marginPx;
            x1 = x2 - scaleBar_px;
            y  = 1 + marginPx;

            textX = (x1 + x2)/2;
            textY = y + textGapPx;

            vAlign = 'top';

        otherwise
            error('Unknown scaleBarCorner option.');
    end

    line(ax, [x1 x2], [y y], 'Color', textColor, 'LineWidth', lineWidth);

    text(ax, textX, textY, sprintf('%g \\mum', scaleBar_um), ...
        'Color', textColor, ...
        'FontSize', fontSize, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', vAlign, ...
        'Interpreter', 'tex');
end

function results = analyzeSignalBackgroundROI(img, cmap, outName, ...
    savePNG, saveFIG, saveMAT, addScaleBar, scaleBar_px, scaleBar_um, ...
    scaleBarCorner, scaleBarMarginPx, scaleBarLineWidth, scaleBarTextColor, ...
    scaleBarFontSize, scaleBarTextGapPx)
% ANALYZESIGNALBACKGROUNDROI  Interactive polygon ROI selection on a RAW
% image to compute SNR and contrast metrics for a reviewer-requested
% quantification.
%
%   results = analyzeSignalBackgroundROI(img, cmap, outName, savePNG, saveFIG, saveMAT, ...
%       addScaleBar, scaleBar_px, scaleBar_um, scaleBarCorner, scaleBarMarginPx, ...
%       scaleBarLineWidth, scaleBarTextColor, scaleBarFontSize, scaleBarTextGapPx)
%
% Draws one polygon ROI for SIGNAL and one for BACKGROUND on img
% (double-click, or right-click > "Finish drawing", to close each
% polygon). Computes:
%   SNR_signal_over_stdSignal      = mean_signal / std_signal
%   SNR_signal_over_stdBackground  = mean_signal / std_background
%   SNR_corrected                  = (mean_signal - mean_background) / std_background
%       (additional, not explicitly requested — standard background-
%        corrected SNR definition in fluorescence microscopy; delete
%        this field/line if you only want the two ratios above)
%   contrast                       = mean_signal / mean_background
%
% Returns a struct with the metrics, the ROI polygon vertex coordinates,
% and the binary masks. If savePNG/saveFIG/saveMAT are true, exports the
% overlay figures (opaque with legend/title, and transparent without) and
% the results struct using outName as the base path.

    figROI = figure('Color','w','Name','Select_SNR_Contrast_ROIs');
    axROI = axes(figROI);
    imagesc(axROI, img, [0 1]);
    axis(axROI, 'image');
    colormap(axROI, cmap);
    colorbar(axROI);

    [nRows, nCols] = size(img);

    title(axROI, 'Draw SIGNAL polygon (double-click to finish) - vertices snap to pixel centers');
    hSignal = drawpolygon(axROI, 'Color', 'r', 'LineWidth', 1.5);
    wait(hSignal);
    snapROIVertices(hSignal, nRows, nCols);   % snap once the polygon is finalized
    addlistener(hSignal, 'ROIMoved', @(src,~) snapROIVertices(src, nRows, nCols));  % keep snapping on later edits
    signalPos  = hSignal.Position;   % [x y] vertex coordinates, snapped to pixel centers
    signalMask = createMask(hSignal, img);
    signalMask = addBoundaryPixels(signalMask, signalPos, nRows, nCols);

    title(axROI, 'Draw BACKGROUND polygon (double-click to finish) - vertices snap to pixel centers');
    hBackground = drawpolygon(axROI, 'Color', 'c', 'LineWidth', 1.5);
    wait(hBackground);
    snapROIVertices(hBackground, nRows, nCols);
    addlistener(hBackground, 'ROIMoved', @(src,~) snapROIVertices(src, nRows, nCols));
    backgroundPos  = hBackground.Position;
    backgroundMask = createMask(hBackground, img);
    backgroundMask = addBoundaryPixels(backgroundMask, backgroundPos, nRows, nCols);

    close(figROI);

    signalVals     = img(signalMask);
    backgroundVals = img(backgroundMask);

    mean_signal     = mean(signalVals);
    std_signal      = std(signalVals);
    mean_background = mean(backgroundVals);
    std_background  = std(backgroundVals);

    SNR_signal_over_stdSignal     = mean_signal / std_signal;
    SNR_signal_over_stdBackground = mean_signal / std_background;
    SNR_corrected = (mean_signal - mean_background) / std_background;
    contrastRatio = mean_signal / mean_background;

    results = struct( ...
        'mean_signal', mean_signal, ...
        'std_signal', std_signal, ...
        'mean_background', mean_background, ...
        'std_background', std_background, ...
        'SNR_signal_over_stdSignal', SNR_signal_over_stdSignal, ...
        'SNR_signal_over_stdBackground', SNR_signal_over_stdBackground, ...
        'SNR_corrected', SNR_corrected, ...
        'contrast', contrastRatio, ...
        'signalROI_xy', signalPos, ...
        'backgroundROI_xy', backgroundPos, ...
        'signalMask', signalMask, ...
        'backgroundMask', backgroundMask);

    % Overlay figure: RAW image with both ROI contours
    figOverlay = figure('Color','w','Name','ROI_Overlay_Raw');
    axOverlay = axes(figOverlay);
    imagesc(axOverlay, img, [0 1]);
    axis(axOverlay, 'image', 'off');
    colormap(axOverlay, cmap);
    colorbar(axOverlay);
    hold(axOverlay, 'on');

    plot(axOverlay, [signalPos(:,1); signalPos(1,1)], ...
        [signalPos(:,2); signalPos(1,2)], 'r-', 'LineWidth', 2);
    plot(axOverlay, [backgroundPos(:,1); backgroundPos(1,1)], ...
        [backgroundPos(:,2); backgroundPos(1,2)], 'c-', 'LineWidth', 2);

    legend(axOverlay, {'Signal ROI','Background ROI'}, ...
        'TextColor', 'k', 'Location', 'southoutside', 'Orientation', 'horizontal');

    title(axOverlay, sprintf( ...
        'SNR_{sig/std_{sig}} = %.2f | SNR_{sig/std_{bg}} = %.2f | Contrast = %.2f', ...
        SNR_signal_over_stdSignal, SNR_signal_over_stdBackground, contrastRatio), ...
        'Interpreter', 'tex');

    if addScaleBar
        drawScaleBar(axOverlay, nRows, nCols, scaleBar_px, scaleBar_um, ...
            scaleBarCorner, scaleBarMarginPx, scaleBarLineWidth, ...
            scaleBarTextColor, scaleBarFontSize, scaleBarTextGapPx);
    end

    drawnow

    if savePNG
        exportgraphics(figOverlay, [outName '_ROI_overlay.png'], 'Resolution', 300);
    end
    if saveFIG
        savefig(figOverlay, [outName '_ROI_overlay.fig']);
    end
    if saveMAT
        save([outName '_SNR_contrast_results.mat'], 'results');
    end

    % ----------------------- FIGURE: RAW + ROI OVERLAY, TRANSPARENT PNG -----------------------
    % Same ROI contours as above, but no scalebar-less version - this one
    % DOES carry the scalebar, just no legend/colorbar/title, transparent
    % background - for direct use in a figure panel.
    figOverlayTransparent = figure('Color','none','Name','ROI_Overlay_Raw_Transparent');
    axOverlayTransparent = axes(figOverlayTransparent);
    imagesc(axOverlayTransparent, img, [0 1]);
    axis(axOverlayTransparent, 'image');
    axis(axOverlayTransparent, 'off');
    colormap(axOverlayTransparent, cmap);
    set(axOverlayTransparent, 'Color', 'none');
    hold(axOverlayTransparent, 'on');

    plot(axOverlayTransparent, [signalPos(:,1); signalPos(1,1)], ...
        [signalPos(:,2); signalPos(1,2)], 'r-', 'LineWidth', 2);
    plot(axOverlayTransparent, [backgroundPos(:,1); backgroundPos(1,1)], ...
        [backgroundPos(:,2); backgroundPos(1,2)], 'c-', 'LineWidth', 2);

    if addScaleBar
        drawScaleBar(axOverlayTransparent, nRows, nCols, scaleBar_px, scaleBar_um, ...
            scaleBarCorner, scaleBarMarginPx, scaleBarLineWidth, ...
            scaleBarTextColor, scaleBarFontSize, scaleBarTextGapPx);
    end

    drawnow

    if savePNG
        exportgraphics(figOverlayTransparent, [outName '_ROI_overlay_transparent.png'], ...
            'BackgroundColor', 'none', 'Resolution', 300);
    end
end

function snapROIVertices(hROI, nRows, nCols)
% SNAPROIVERTICES  Snap all vertices of a drawpolygon ROI to the nearest
% pixel center. images.roi.Polygon has no event that fires per-vertex
% during initial interactive drawing, so this is called once explicitly
% right after wait(hROI) returns (i.e. once the polygon is finalized by
% double-click), and again via the 'ROIMoved' listener if the user drags
% a vertex afterward. Relies on imagesc's default convention: pixel
% (row=i, col=j) is centered at axes coordinates (x=j, y=i), so snapping
% is a simple round.

    pos = hROI.Position;
    pos(:,1) = min(max(round(pos(:,1)), 1), nCols);   % x -> column index
    pos(:,2) = min(max(round(pos(:,2)), 1), nRows);   % y -> row index
    hROI.Position = pos;
end

function mask = addBoundaryPixels(mask, polyXY, nRows, nCols)
% ADDBOUNDARYPIXELS  OR the pixels lying exactly on the polygon boundary
% into an existing interior mask. createMask/poly2mask can inconsistently
% include/exclude pixels sitting exactly on an edge (a scan-fill tie-break
% issue that becomes common once vertices are snapped to pixel centers,
% since edges then frequently run exactly along a row/column of centers).
% Each edge is rasterized with Bresenham's algorithm so every pixel the
% drawn boundary line passes through is guaranteed to be included.

    nV = size(polyXY, 1);
    for i = 1:nV
        p1 = polyXY(i, :);
        p2 = polyXY(mod(i, nV) + 1, :);   % wrap last vertex back to first
        pts = bresenhamLine(p1(1), p1(2), p2(1), p2(2));
        pts(:,1) = min(max(pts(:,1), 1), nCols);
        pts(:,2) = min(max(pts(:,2), 1), nRows);
        idx = sub2ind([nRows, nCols], pts(:,2), pts(:,1));
        mask(idx) = true;
    end
end

function pts = bresenhamLine(x1, y1, x2, y2)
% BRESENHAMLINE  Integer pixel coordinates along the line from
% (x1,y1) to (x2,y2), using Bresenham's line algorithm. Endpoints are
% assumed to already be integer pixel indices (as guaranteed by
% snapROIVertices).

    x1 = round(x1); y1 = round(y1);
    x2 = round(x2); y2 = round(y2);

    dx = abs(x2 - x1); dy = abs(y2 - y1);
    sx = sign(x2 - x1); sy = sign(y2 - y1);
    err = dx - dy;

    x = x1; y = y1;
    pts = zeros(0, 2);
    while true
        pts(end+1, :) = [x, y]; %#ok<AGROW>
        if x == x2 && y == y2
            break
        end
        e2 = 2*err;
        if e2 > -dy
            err = err - dy;
            x = x + sx;
        end
        if e2 < dx
            err = err + dx;
            y = y + sy;
        end
    end
end