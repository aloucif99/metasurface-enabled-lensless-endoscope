% Mode-overlap efficiency versus rotational MS-MCF misalignment.

clear
close all
clc

%% Parameters

lambda      = 920e-9;      % [m]
period_MS  = 0.4e-6;       % [m]
dx_ms      = period_MS;    % [m]
n_cores    = 120;

M_mag          = 1;
size_MS_x_um   = 180;
size_MS_y_um   = 180;

Nx_MS = round((size_MS_x_um*1e-6) / period_MS);
Ny_MS = Nx_MS;
N_MS  = Nx_MS;

% Fixed z0 and propagation distance
z0_MS  = -30e-6;     % [m]
f_fixed = 28e-6;     % [m]

% Fixed lateral offset
x_offset = 0;        % [m]
y_offset = 0;        % [m]

% Rotation sweep
rotation_deg_list = 0:0.5:5;
num_angles = numel(rotation_deg_list);

% Rotation angles shown as core coupling maps
rotation_display_deg = [0 1 2 5];
num_display = numel(rotation_display_deg);

% MS-plane coordinates
x_ms = dx_ms*((-N_MS/2):(N_MS/2-1));
y_ms = x_ms;

[X_ms, Y_ms] = meshgrid(x_ms, y_ms);

%% Phase vector

load(fullfile(fileparts(mfilename('fullpath')), 'PhaseVector400um.mat'),'phase_vec');
Phi = phase_vec(:);
Phi = Phi(1:n_cores);

%% Fermat core positions

p_MS = 30;
size_beamlet_ms = p_MS / 1.72;

rho      = sqrt(1:1000)';
theta_sp = (3 - sqrt(5)) * pi * (1:1000)';

poscore_ms = size_beamlet_ms * ...
    [rho .* cos(theta_sp), rho .* sin(theta_sp)];

pitch_ms    = p_MS * period_MS;
pitch_fiber = pitch_ms / M_mag;

%% Phase mask input

Amp      = ones(n_cores,1);
basevecX = Amp .* exp(1i * Phi);

rot_angle = 0;
M         = 1;
theta_rot = 0;
modPx     = Inf;
modPy     = Inf;
offsetx   = 0;
offsety   = 0;
twopigray = 256;
save2tif  = [];
output2WS = true;
basedir   = [];
basefname = [];

%% Fiber grid and mode

dx_fiber = 400e-9;                        % [m]
size_fiber_x = size_MS_x_um * 1e-6;       % [m]
size_fiber_y = size_MS_y_um * 1e-6;       % [m]

Nx_fiber = round(size_fiber_x / dx_fiber);
Ny_fiber = round(size_fiber_y / dx_fiber);
N_fiber  = Nx_fiber;

x_fiber = dx_fiber*((-N_fiber/2):(N_fiber/2-1));
y_fiber = x_fiber;

[X_fiber, Y_fiber] = meshgrid(x_fiber, y_fiber);

%% Fiber parameters and Gaussian mode using GRIN model

CoreDiam = 2.5e-6;       % [m] core diameter
a        = CoreDiam/2;   % [m] core radius

ncoeur = 1.48;
ngaine = 1.45;
g_grin = 2;              % Parabolic GRIN profile

k0 = 2*pi/lambda;

NA_core = sqrt(ncoeur^2 - ngaine^2);

Delta_grin = (ncoeur^2 - ngaine^2) / (2*ncoeur^2);

V = k0 * ncoeur * a * sqrt(2*Delta_grin);

A_grin = sqrt(2/5 * (1 + 4*(2/g_grin)^(5/6)));
B_grin = exp(0.298/g_grin) - 1 ...
       + 1.478*(1 - exp(-0.077*g_grin));
C_grin = 3.76 + exp(4.19/g_grin^0.418);

wg_fiber = a * ( ...
    A_grin / V^(2/(g_grin+2)) + ...
    B_grin / V^(3/2) + ...
    C_grin / V^6 );

MFD_fiber = 2 * wg_fiber;

NA_mode_estimated = sin(atan(lambda/(pi*wg_fiber)));
%% Gaussian fiber mode on fiber grid

F_fiber = exp(-(X_fiber.^2 + Y_fiber.^2)/wg_fiber^2);
F_fiber = F_fiber ./ sqrt(sum(sum(abs(F_fiber).^2)));

F_fft = fftshift(fft2(ifftshift(F_fiber)));

%% Core positions in fiber plane + segmentation

poscore_fiber_phys = poscore_ms * period_MS;   % [m]
poscore_fiber_pix  = poscore_fiber_phys / dx_fiber;

[IN_fiber, SegSizes_fiber] = local_voronoi_segmentation( ...
    poscore_fiber_pix, Nx_fiber, Ny_fiber);

NpixTot_fiber = N_fiber * N_fiber;

%% PART 5 — FOURIER GRID ON MS PLANE

fx_ms = 1/(N_MS*dx_ms)*((-N_MS/2):(N_MS/2-1));
fy_ms = fx_ms;

[Fx_ms, FY_ms] = meshgrid(fx_ms, fy_ms);

%% PART 6 — GENERATE FIXED MS PHASE MASK

phi2_MS = pi/(lambda*z0_MS) * dx_ms^2;

test_ms_slm = phase_mask_mcf(basevecX, poscore_ms, phi2_MS, Nx_MS, Ny_MS, rot_angle);

test_ms_slm = squeeze(test_ms_slm);

%% PART 7 — OPTIONAL DISPLAY OF PHASE MASK

figure('Color','w','Position',[200 200 650 550]);

imagesc(x_ms*1e6, y_ms*1e6, test_ms_slm);
axis image;

colormap(jet);
colorbar;

xlabel('$x_{\mathrm{MS}}~(\mu\mathrm{m})$', 'Interpreter','latex');
ylabel('$y_{\mathrm{MS}}~(\mu\mathrm{m})$', 'Interpreter','latex');

title(sprintf('MS phase mask, $z_0 = %.0f~\\mu\\mathrm{m}$', abs(z0_MS)*1e6), ...
    'Interpreter','latex', ...
    'FontWeight','normal');

set(gca, ...
    'FontSize',16, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.2);

%% PART 8 — PROPAGATE FIELD TO FIXED DISTANCE

tempfield_ms = exp(1i * test_ms_slm);
tempfield_fft_ms = fftshift(fft2(ifftshift(tempfield_ms)));

H_ms = exp(-1i * pi * lambda * f_fixed .* (Fx_ms.^2 + FY_ms.^2));

spotfield_fft_ms = tempfield_fft_ms .* H_ms;
spotfield_ms = fftshift(ifft2(ifftshift(spotfield_fft_ms)));

%% PART 9 — OPTIONAL DISPLAY OF INTENSITY AT FIXED DISTANCE

I_ms = abs(spotfield_ms).^2;
I_ms = I_ms ./ max(I_ms(:));

figure('Color','w','Position',[250 250 650 550]);

imagesc(x_ms*1e6, y_ms*1e6, I_ms);
axis image;

colormap(turbo);
colorbar;

xlabel('$x_{\mathrm{MS}}~(\mu\mathrm{m})$', 'Interpreter','latex');
ylabel('$y_{\mathrm{MS}}~(\mu\mathrm{m})$', 'Interpreter','latex');

title(sprintf('Spot pattern at $z = %.0f~\\mu\\mathrm{m}$', f_fixed*1e6), ...
    'Interpreter','latex', ...
    'FontWeight','normal');

set(gca, ...
    'FontSize',16, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.2);

%% PART 10 — SWEEP ROTATION ANGLE

results_overlap_all = NaN(n_cores, num_angles);

mean_overlap_vs_rot = NaN(num_angles, 1);
std_overlap_vs_rot  = NaN(num_angles, 1);

coupling_maps_display = NaN(n_cores, num_display);

for ai = 1:num_angles

    rot_deg = rotation_deg_list(ai);
    rot_rad = rot_deg * pi / 180;

    % ---------------------------------------------------------------------
    % Apply rotational mismatch between MS-generated spot pattern and
    % fixed MCF core positions.
    %
    % x_offset = 0
    % y_offset = 0
    % ---------------------------------------------------------------------

    Xc = X_fiber - x_offset;
    Yc = Y_fiber - y_offset;

    Xq = M_mag * ( cos(rot_rad)*Xc + sin(rot_rad)*Yc );
    Yq = M_mag * (-sin(rot_rad)*Xc + cos(rot_rad)*Yc );

    spotfield_fiber = interp2( ...
        X_ms, Y_ms, spotfield_ms, Xq, Yq, 'linear', 0);

    % Normalize field on fiber grid
    norm_spot = sqrt(sum(sum(abs(spotfield_fiber).^2)));

    if norm_spot > 0
        spotfield_fiber = spotfield_fiber / norm_spot;
    end

    % ---------------------------------------------------------------------
    % Overlap using convolution theorem
    % ---------------------------------------------------------------------

    S_fft_fiber = fftshift(fft2(ifftshift(spotfield_fiber)));

    C_fft = F_fft .* S_fft_fiber;
    C = fftshift(ifft2(ifftshift(C_fft)));

    % ---------------------------------------------------------------------
    % Coupling efficiency:
    % evaluate at fixed MCF core centers
    % ---------------------------------------------------------------------

    coupleff = NaN(n_cores, 1);

    for k = 1:n_cores

        x_core = poscore_fiber_phys(k,1);
        y_core = poscore_fiber_phys(k,2);

        [~, col_core] = min(abs(x_fiber - x_core));
        [~, row_core] = min(abs(y_fiber - y_core));

        if isempty(IN_fiber{k}) || SegSizes_fiber(k) == 0
            coupleff(k) = NaN;
            continue;
        end

        coupleff(k) = abs(C(row_core, col_core)).^2 * ...
            (NpixTot_fiber / SegSizes_fiber(k));

    end

    results_overlap_all(:, ai) = coupleff * 100;

    mean_overlap_vs_rot(ai) = mean(results_overlap_all(:,ai), 'omitnan');
    std_overlap_vs_rot(ai)  = std(results_overlap_all(:,ai), [], 'omitnan');

    % Store selected angles for core-map display
    idx_disp = find(abs(rotation_display_deg - rot_deg) < 1e-12);

    if ~isempty(idx_disp)
        coupling_maps_display(:, idx_disp) = results_overlap_all(:, ai);
    end

end

%% PART 11 — FIGURE 1: MEAN OVERLAP VS ROTATION

fig1 = figure('Color','w','Position',[120 120 980 530]);

errorbar(rotation_deg_list, mean_overlap_vs_rot, std_overlap_vs_rot, ...
    '-o', ...
    'LineWidth',2, ...
    'MarkerSize',7, ...
    'CapSize',7);

grid on
box on

xlabel('Rotation angle $(^\circ)$', ...
    'Interpreter','latex', ...
    'FontSize',20);

ylabel('Mean overlap efficiency (\%)', ...
    'Interpreter','latex', ...
    'FontSize',20);

title(sprintf('Overlap efficiency vs rotational MS - MCF misalignment', ...
    f_fixed*1e6), ...
    'Interpreter','latex', ...
    'FontSize',23, ...
    'FontWeight','normal');

legend('Mean $\pm$ std over 120 cores', ...
    'Interpreter','latex', ...
    'Location','best', ...
    'FontSize',20);

xlim([0 5]);
xticks(0:0.5:5);
ylim([0 100]);

set(gca, ...
    'FontSize',22, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.2, ...
    'TickDir','out');

%% PART 12 — FIGURE 2: CORE MAPS COLORED BY OVERLAP EFFICIENCY

fig2 = figure('Color','w','Position',[80 80 800 800]);

tiledlayout(2,2, ...
    'Padding','compact', ...
    'TileSpacing','compact');

% Use only the first 120 cores because overlap is calculated for n_cores only
x_core_um = poscore_fiber_phys(1:n_cores,1) * 1e6;
y_core_um = poscore_fiber_phys(1:n_cores,2) * 1e6;

for di = 1:num_display

    ax = nexttile;

    colorValues = coupling_maps_display(:,di);
    colorValues = colorValues(:);

    scatter(ax, x_core_um, y_core_um, 95, colorValues, 'filled');

    axis(ax, 'image');

    xlim(ax, [-90 90]);
    ylim(ax, [-90 90]);

    colormap(ax, 'parula');
    caxis(ax, [0 100]);

    xlabel(ax, '$x~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',20);

    ylabel(ax, '$y~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',20);

    title(ax, sprintf('$\\theta = %g^{\\circ}$', rotation_display_deg(di)), ...
        'Interpreter','latex', ...
        'FontSize',23, ...
        'FontWeight','normal');

    set(ax, ...
        'FontSize',20, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1.1, ...
        'TickDir','out', ...
        'Box','on');

end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 20;
cb.Label.String = 'Mode overlap efficiency (\%)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 22;
