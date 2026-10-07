% Distal-facet core intensity fraction from the supplied MS and SLM .fig files.
% All paths are relative to this script. Requires Image Processing Toolbox.

clear; close all; clc
folder = fileparts(mfilename('fullpath'));
outputDir = fullfile(folder, 'output');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end

caseNames = {'MS', 'SLM'};
figNames = {'DistalEnd_MS.fig', 'DistalEnd_SLM.fig'};
NcoresExpected = 120;
% [x center, y center, radius] in image pixels; NaN = select 3 boundary points.
% The original MS selection was not saved, so do not invent its coordinates.
fiberCircleByCase = [NaN NaN NaN; 590 615 430];

% Detection and integration parameters from the original analysis.
minPeakDist_px = 10;
coreSmoothSigma_px = 1.2;
backgroundSigma_px = 25;
peakMinRel = 0.006;
peakSigmaThreshold = 1.5;
maxThresholdRel = 0.020;
fwhmSearchRadius_px = 14;
minFWHM_px = 2;
maxFWHM_px = 22;
fwhmApertureScale = 2; % aperture diameter = 2 x equivalent FWHM diameter

results = struct([]);
imagesForDisplay = cell(1,2);
for icase = 1:numel(caseNames)
    name = caseNames{icase};
    sourceFile = fullfile(folder, figNames{icase});
    assert(isfile(sourceFile), 'Missing input FIG file: %s', sourceFile);
    h = openfig(sourceFile, 'invisible');
    imageHandles = findobj(h, 'Type', 'image');
    assert(~isempty(imageHandles), 'No image object found in %s', sourceFile);
    sizes = arrayfun(@(obj) numel(get(obj,'CData')), imageHandles);
    [~, largest] = max(sizes);
    Iraw = double(get(imageHandles(largest), 'CData'));
    close(h);
    assert(ismatrix(Iraw), '%s must contain a two-dimensional intensity image.', sourceFile);
    Iraw(~isfinite(Iraw)) = 0;
    [Ny, Nx] = size(Iraw);

    % Fixed detector offset determined from the image border.
    border = max(1, round(0.05*min(Ny,Nx)));
    bgPixels = [reshape(Iraw(1:border,:),[],1); ...
                reshape(Iraw(end-border+1:end,:),[],1); ...
                reshape(Iraw(:,1:border),[],1); ...
                reshape(Iraw(:,end-border+1:end),[],1)];
    backgroundOffset = median(bgPixels);
    I = max(Iraw-backgroundOffset, 0);
    [X,Y] = meshgrid(1:Nx,1:Ny);

    circle = fiberCircleByCase(icase,:);
    if ~all(isfinite(circle))
        fig = figure('Color','w'); imagesc(I); axis image; colormap(parula);
        title(sprintf('%s: select 3 points along the real fiber boundary', name));
        xlabel('x pixel'); ylabel('y pixel');
        [cx,cy] = ginput(3); close(fig);
        assert(numel(cx)==3, 'Three boundary points are required.');
        coeff = [2*cx(:), 2*cy(:), ones(3,1)] \ (cx(:).^2+cy(:).^2);
        circle = [coeff(1),coeff(2),sqrt(coeff(3)+coeff(1)^2+coeff(2)^2)];
        fprintf('%s fiber boundary, for reproducible reruns: [%.3f %.3f %.3f]\n', name,circle);
        fprintf('Set fiberCircleByCase(%d,:) to this value after inspecting the overlay.\n',icase);
    end
    xc=circle(1); yc=circle(2); R=circle(3);
    assert(isfinite(R) && R>0, 'Invalid fiber radius for %s.',name);
    fiberMask = (X-xc).^2 + (Y-yc).^2 <= R^2;

    % Find bright local maxima in a high-pass copy; integrate ORIGINAL
    % background-corrected intensities (not the high-pass detection image).
    Ismooth = imgaussfilt(I, coreSmoothSigma_px);
    Ihigh = max(Ismooth-imgaussfilt(Ismooth,backgroundSigma_px),0);
    Ihigh(~fiberMask)=0;
    Idetect=mat2gray(Ihigh);
    assert(any(Idetect(:)>0), 'Detection image is empty: check %s fiber boundary.',name);
    values=Idetect(fiberMask); values=values(isfinite(values));
    medVal=median(values);
    sigmaRobust=1.4826*median(abs(values-medVal));
    if sigmaRobust==0, sigmaRobust=std(values); end
    thr = min(max(peakMinRel,medVal+peakSigmaThreshold*sigmaRobust),maxThresholdRel);
    if thr<=0 || thr>=1, thr=peakMinRel; end
    localMax=imregionalmax(Idetect) & fiberMask & Idetect>=thr;
    [yy,xx]=find(localMax);
    strengths=Idetect(sub2ind(size(Idetect),yy,xx));
    [~,order]=sort(strengths,'descend');
    xx=xx(order); yy=yy(order);
    candidateCenters=zeros(0,2);
    for k=1:numel(xx)
        pt=[xx(k),yy(k)];
        if isempty(candidateCenters) || all(vecnorm(candidateCenters-pt,2,2)>=minPeakDist_px)
            candidateCenters(end+1,:)=pt; %#ok<SAGROW>
        end
    end

    % The optional file records only additional, visually CONFIRMED cores
    % missed by the automatic detector; columns are x_px,y_px. This preserves
    % an auditable manual correction without forcing artificial detections.
    reviewedFile=fullfile(folder,sprintf('reviewed_extra_centers_%s.csv',name));
    additional=zeros(0,2);
    if isfile(reviewedFile)
        additional=readmatrix(reviewedFile);
        assert(size(additional,2)==2 && all(isfinite(additional(:))), ...
            '%s must have exactly two numeric columns: x_px,y_px.',reviewedFile);
    end
    allCandidates=[candidateCenters; additional];
    isManual=[false(size(candidateCenters,1),1);true(size(additional,1),1)];
    centers=zeros(0,2); coreFWHM_px=zeros(0,1); acceptedManual=false(0,1);
    coreMask=false(Ny,Nx);
    rejected=zeros(0,2);
    for k=1:size(allCandidates,1)
        x0=round(allCandidates(k,1)); y0=round(allCandidates(k,2));
        if x0<1 || x0>Nx || y0<1 || y0>Ny || ~fiberMask(y0,x0)
            rejected(end+1,:)=[x0 y0]; %#ok<SAGROW>
            if isManual(k), warning('%s reviewed point %d is outside the fiber.',name,k); end
            continue;
        end
        if ~isempty(centers) && any(vecnorm(centers-[x0 y0],2,2)<minPeakDist_px)
            if isManual(k), warning('%s reviewed point [%.0f %.0f] duplicates a core.',name,x0,y0); end
            continue;
        end
        x1=max(1,x0-fwhmSearchRadius_px); x2=min(Nx,x0+fwhmSearchRadius_px);
        y1=max(1,y0-fwhmSearchRadius_px); y2=min(Ny,y0+fwhmSearchRadius_px);
        patch=Ismooth(y1:y2,x1:x2);
        localSignal=max(patch-median(patch(:)),0);
        localX=x0-x1+1; localY=y0-y1+1;
        peak=localSignal(localY,localX);
        diam=NaN;
        if peak>0
            components=bwconncomp(localSignal>=0.5*peak);
            pixelIndex=sub2ind(size(localSignal),localY,localX);
            for j=1:components.NumObjects
                if any(components.PixelIdxList{j}==pixelIndex)
                    diam=2*sqrt(numel(components.PixelIdxList{j})/pi);
                    break;
                end
            end
        end
        if ~isfinite(diam) || diam<minFWHM_px || diam>maxFWHM_px
            rejected(end+1,:)=[x0 y0]; %#ok<SAGROW>
            if isManual(k)
                warning('%s reviewed core [%.0f %.0f] has invalid FWHM; inspect its center.',name,x0,y0);
            end
            continue;
        end
        radius=fwhmApertureScale*diam/2;
        coreMask=coreMask | (fiberMask & ((X-x0).^2+(Y-y0).^2<=radius^2));
        centers(end+1,:)=[x0 y0]; %#ok<SAGROW>
        coreFWHM_px(end+1,1)=diam; %#ok<SAGROW>
        acceptedManual(end+1,1)=isManual(k); %#ok<SAGROW>
        if size(centers,1)>=NcoresExpected, break; end
    end

    nDetected=size(centers,1);
    if nDetected~=NcoresExpected
        warning('%s: %d / %d accepted. Inspect the saved overlay; do not assume 120 detections.', ...
            name,nDetected,NcoresExpected);
    end
    sumCore=sum(I(coreMask)); sumFacet=sum(I(fiberMask));
    assert(sumFacet>0,'%s: whole-facet intensity is zero.',name);
    results(icase).configuration=name;
    results(icase).nominalPhysicalCoreCount=NcoresExpected;
    results(icase).acceptedCoreCount=nDetected;
    results(icase).candidateCount=size(candidateCenters,1);
    results(icase).reviewedExtraCount=nnz(acceptedManual);
    results(icase).fiberCircle=circle;
    results(icase).backgroundOffset=backgroundOffset;
    results(icase).threshold=thr;
    results(icase).centers=centers;
    results(icase).coreFWHM_px=coreFWHM_px;
    results(icase).fiberMask=fiberMask;
    results(icase).coreMask=coreMask;
    results(icase).coreIntensity=sumCore;
    results(icase).facetIntensity=sumFacet;
    results(icase).coreFraction=sumCore/sumFacet;
    results(icase).claddingFraction=1-results(icase).coreFraction;
    imagesForDisplay{icase}=I;
    writematrix(centers,fullfile(outputDir,sprintf('detected_core_centers_%s.csv',name)));
    fprintf('%s: %d thresholded candidates, %d accepted cores, %.3f%% core fraction.\n', ...
        name,size(candidateCenters,1),nDetected,100*results(icase).coreFraction);
end

T=table({results.configuration}',[results.nominalPhysicalCoreCount]', ...
    [results.acceptedCoreCount]',100*[results.coreFraction]', ...
    100*[results.claddingFraction]', ...
    'VariableNames',{'Configuration','PhysicalCoreCount','AcceptedCores', ...
    'CoreFraction_percent','OutsideCoreFraction_percent'});
disp(T);
save(fullfile(outputDir,'core_fraction_results.mat'),'results','T');
f=figure('Color','w'); tiledlayout(1,2,'TileSpacing','compact');
for icase=1:2
    nexttile; imagesc(imagesForDisplay{icase}); axis image; colormap(parula); hold on;
    contour(double(results(icase).fiberMask),[0.5 0.5],'w','LineWidth',1);
    contour(double(results(icase).coreMask),[0.5 0.5],'m','LineWidth',0.5);
    plot(results(icase).centers(:,1),results(icase).centers(:,2), ...
        'c.','MarkerSize',5);
    title(sprintf('%s: %.2f%% (%d of %d cores)',results(icase).configuration, ...
        100*results(icase).coreFraction,results(icase).acceptedCoreCount,NcoresExpected));
    xlabel('x (pixel)'); ylabel('y (pixel)');
end
exportgraphics(f,fullfile(outputDir,'core_fraction_diagnostic.png'),'Resolution',200);
