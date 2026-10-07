% Z-scan analysis of the metasurface focal array.


clear; clc; close all;

%% ---------------- PARAMETERS ----------------

dz_um = 0.2;                 % z step [um]
z_ref_um = 30;               % reference plane [um]

pixel_size_um = 3.45/20;     % camera pixel / magnification [um]
aperture_um = 12;            % width of X-Z cut [um]

nLensesToPlot = 10;

threshold = 0.30;            % detection threshold
kernel_size = 10;            % local integration for detection

avg_half_height = 2;         % average over +/- 2 rows
interpFactor = 6;            % interpolation along z
z_plot_max_um = 60;


%% ---------------- LOAD DATA ----------------

frames = double(h5read(fullfile(fileparts(mfilename('fullpath')), 'data.hdf5'),'/frames'));

[Ny,Nx,Nz] = size(frames);

z_um = (0:Nz-1)*dz_um;


%% ---------------- REFERENCE PLANE ----------------

% Acquired plane closest to z_ref_um
[~,idxRef] = min(abs(z_um-z_ref_um));

z_ref_actual_um = z_um(idxRef);

ref = frames(:,:,idxRef);


%% ---------------- DETECT ALL LENSLETS ----------------

% Local integration used only for robust spot detection
kernel = ones(kernel_size,kernel_size);

ref_detect = conv2(ref,kernel,'same');


% Detect local maxima
BW = imregionalmax(ref_detect);

BW = BW & ...
    (ref_detect > threshold*max(ref_detect(:)));


% Remove image borders
margin = 20;

BW(1:margin,:) = 0;
BW(end-margin+1:end,:) = 0;

BW(:,1:margin) = 0;
BW(:,end-margin+1:end) = 0;


% Coordinates of all detected lenslets
[rowPeak,colPeak] = find(BW);

nPeaks = length(rowPeak);

fprintf('Detected lenslets: %d\n',nPeaks);


%% ---------------- SHOW REFERENCE IMAGE ----------------

figure;

imagesc(ref);

axis image;
colormap hot;
colorbar;

hold on;

% Show all detected lenslets
plot( ...
    colPeak,rowPeak, ...
    'co', ...
    'MarkerSize',8, ...
    'LineWidth',1.2);

title(sprintf( ...
    'Reference plane: z = %.1f \\mum, N = %d', ...
    z_ref_actual_um,nPeaks));


%% ---------------- X-Z CUTS OF FIRST 10 LENSLETS ----------------

nLensesToPlot = min(nLensesToPlot,nPeaks);

half_width_px = round( ...
    (aperture_um/2)/pixel_size_um );

z_fine = linspace( ...
    z_um(1), ...
    z_um(end), ...
    Nz*interpFactor );


for p = 1:nLensesToPlot

    % First detected lenslets
    r0 = rowPeak(p);
    c0 = colPeak(p);


    %% X window around lenslet

    c1 = max(1,c0-half_width_px);
    c2 = min(Nx,c0+half_width_px);

    x_um = ((c1:c2)-c0)*pixel_size_um;


    %% Build X-Z map

    xz_map = zeros(Nz,c2-c1+1);


    for iz = 1:Nz

        r1 = max(1,r0-avg_half_height);
        r2 = min(Ny,r0+avg_half_height);

        patch = frames(r1:r2,c1:c2,iz);

        % Average over a few pixels in Y
        xz_map(iz,:) = mean(patch,1);

    end


    %% Normalize

    if max(xz_map(:)) > 0
        xz_map = xz_map/max(xz_map(:));
    end


    %% Interpolate along Z

    xz_map_fine = interp1( ...
        z_um, ...
        xz_map, ...
        z_fine, ...
        'pchip');


    %% Plot

    figure('Position',[100 100 300 480]);

    imagesc( ...
        x_um, ...
        z_fine, ...
        xz_map_fine);

    axis xy;

    colormap hot;
    colorbar;

    caxis([0 1]);

    ylim([0 min(z_plot_max_um,z_um(end))]);


    xlabel('$x~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',24);

    ylabel('$z~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',24);


    set(gca, ...
        'FontSize',20, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1.2);

end