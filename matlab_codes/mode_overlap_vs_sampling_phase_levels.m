% Mode-overlap efficiency versus lateral sampling pitch and phase discretization.
% Requires PhaseVector400um.mat in the same folder.

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

z0_MS  = -30e-6;     % [m]
f_fixed = 28e-6;     % [m]

x_offset = 0;
y_offset = 0;

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

%% Fiber grid and mode

dx_fiber = 400e-9;

size_fiber_x = size_MS_x_um * 1e-6;
size_fiber_y = size_MS_y_um * 1e-6;

Nx_fiber = round(size_fiber_x / dx_fiber);
Ny_fiber = Nx_fiber;
N_fiber  = Nx_fiber;

x_fiber = dx_fiber*((-N_fiber/2):(N_fiber/2-1));
y_fiber = x_fiber;

[X_fiber, Y_fiber] = meshgrid(x_fiber, y_fiber);

CoreDiam = 2.5e-6;
a        = CoreDiam/2;

ncoeur = 1.48;
ngaine = 1.45;

g_grin = 2;

k0 = 2*pi/lambda;

Delta_grin = (ncoeur^2 - ngaine^2) / (2*ncoeur^2);
V = k0 * ncoeur * a * sqrt(2*Delta_grin);

A_grin = sqrt(2/5 * (1 + 4 * (2/g_grin)^(5/6)));
B_grin = exp(0.298/g_grin) - 1 + 1.478 * (1 - exp(-0.077*g_grin));
C_grin = 3.76 + exp(4.19 / g_grin^(0.418));

wg_fiber = a * ( ...
    A_grin / V^(2/(g_grin+2)) + ...
    B_grin / V^(3/2) + ...
    C_grin / V^6 );

MFD_fiber = 2 * wg_fiber;

F_fiber = exp(-(X_fiber.^2 + Y_fiber.^2)/wg_fiber^2);
F_fiber = F_fiber ./ sqrt(sum(sum(abs(F_fiber).^2)));

F_fft = fftshift(fft2(ifftshift(F_fiber)));

poscore_fiber_phys = poscore_ms * period_MS;
poscore_fiber_pix  = poscore_fiber_phys / dx_fiber;

[IN_fiber, SegSizes_fiber] = local_voronoi_segmentation( ...
    poscore_fiber_pix, Nx_fiber, Ny_fiber);

NpixTot_fiber = N_fiber * N_fiber;

%% Fourier grid

fx_ms = 1/(N_MS*dx_ms)*((-N_MS/2):(N_MS/2-1));
fy_ms = fx_ms;

[Fx_ms, FY_ms] = meshgrid(fx_ms, fy_ms);

%% Continuous phase mask

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

phi2_MS = pi/(lambda*z0_MS) * dx_ms^2;

phase_continuous = phase_mask_mcf(basevecX, poscore_ms, phi2_MS, Nx_MS, Ny_MS, rot_angle);

phase_continuous = squeeze(phase_continuous);
phase_continuous = mod(phase_continuous, 2*pi);

%% PART 6 — CASES

samplePitch_um_list = [3.0 2.0 1.0 0.4];

samplePitchLabels = { ...
    '3 $\mu$m sampling', ...
    '2 $\mu$m sampling', ...
    '1 $\mu$m sampling', ...
    '400 nm sampling'};

caseNames = {'2 levels', '3 levels', '6 levels', '9 levels', '12 levels', 'Continuous'};
Nlevels_list = [2, 3, 6, 9, 12, Inf];

n_pitch = numel(samplePitch_um_list);
n_cases = numel(Nlevels_list);

plotCaseNames = {'3 levels', '12 levels', 'Continuous'};
plotCaseIdx = [2 5 6];
n_plot_cases = numel(plotCaseIdx);

phase_masks = zeros(Ny_MS, Nx_MS, n_pitch, n_cases);

mean_overlap = NaN(n_pitch, n_cases);
std_overlap  = NaN(n_pitch, n_cases);

I_fiber_cases = zeros(N_fiber, N_fiber, n_pitch, n_cases, 'single');

%% PART 7 — LOOP OVER SAMPLING PITCH AND PHASE LEVELS

for pi_idx = 1:n_pitch

    sample_pitch_um = samplePitch_um_list(pi_idx);
    sample_pitch = sample_pitch_um * 1e-6;

    if sample_pitch <= dx_ms * 1.01

        phase_sampled = phase_continuous;

    else

        X_sample = sample_pitch * round(X_ms / sample_pitch);
        Y_sample = sample_pitch * round(Y_ms / sample_pitch);

        phase_sampled = interp2( ...
            X_ms, Y_ms, phase_continuous, ...
            X_sample, Y_sample, ...
            'nearest', 0);

        phase_sampled = mod(phase_sampled, 2*pi);

    end

    for ci = 1:n_cases

        if isinf(Nlevels_list(ci))

            phase_case = phase_sampled;

        else

            Nlevels = Nlevels_list(ci);
            dphi = 2*pi / Nlevels;

            level_index = round(phase_sampled / dphi);
            level_index = mod(level_index, Nlevels);

            phase_case = level_index * dphi;

        end

        phase_case = mod(phase_case, 2*pi);
        phase_masks(:,:,pi_idx,ci) = phase_case;

        tempfield_ms = exp(1i * phase_case);

        tempfield_fft_ms = fftshift(fft2(ifftshift(tempfield_ms)));

        H_ms = exp(-1i * pi * lambda * f_fixed .* ...
            (Fx_ms.^2 + FY_ms.^2));

        spotfield_fft_ms = tempfield_fft_ms .* H_ms;
        spotfield_ms = fftshift(ifft2(ifftshift(spotfield_fft_ms)));

        Xq = M_mag * (X_fiber - x_offset);
        Yq = M_mag * (Y_fiber - y_offset);

        spotfield_fiber = interp2( ...
            X_ms, Y_ms, spotfield_ms, Xq, Yq, 'linear', 0);

        norm_spot = sqrt(sum(sum(abs(spotfield_fiber).^2)));

        if norm_spot > 0
            spotfield_fiber = spotfield_fiber / norm_spot;
        end

        I_fiber = abs(spotfield_fiber).^2;

        if max(I_fiber(:)) > 0
            I_fiber_cases(:,:,pi_idx,ci) = single(I_fiber ./ max(I_fiber(:)));
        else
            I_fiber_cases(:,:,pi_idx,ci) = single(I_fiber);
        end

        S_fft_fiber = fftshift(fft2(ifftshift(spotfield_fiber)));

        C_fft = F_fft .* S_fft_fiber;
        C = fftshift(ifft2(ifftshift(C_fft)));

        mode_overlap = NaN(n_cores,1);

        for k = 1:n_cores

            x_core = poscore_fiber_phys(k,1);
            y_core = poscore_fiber_phys(k,2);

            [~, col_core] = min(abs(x_fiber - x_core));
            [~, row_core] = min(abs(y_fiber - y_core));

            if isempty(IN_fiber{k}) || SegSizes_fiber(k) == 0
                mode_overlap(k) = NaN;
                continue;
            end

            mode_overlap(k) = abs(C(row_core, col_core)).^2 * ...
                (NpixTot_fiber / SegSizes_fiber(k));

        end

        mode_overlap = mode_overlap * 100;

        mean_overlap(pi_idx,ci) = mean(mode_overlap, 'omitnan');
        std_overlap(pi_idx,ci)  = std(mode_overlap, [], 'omitnan');

    end

end

%% FIGURE 1 — PHASE MASKS BEFORE PROPAGATION

fig1 = figure('Color','w', 'Units','pixels', 'Position',[40 60 1350 1000]);

t1 = tiledlayout(n_pitch, n_plot_cases, ...
    'Padding','compact', ...
    'TileSpacing','compact');

for pi_idx = 1:n_pitch

    for pci = 1:n_plot_cases

        ci = plotCaseIdx(pci);

        ax = nexttile;

        imagesc(ax, x_ms*1e6, y_ms*1e6, phase_masks(:,:,pi_idx,ci));

        axis(ax, 'image');

        xlim(ax, [-90 90]);
        ylim(ax, [-90 90]);

        set(ax, 'YDir', 'normal');

        colormap(ax, turbo);
        caxis(ax, [0 2*pi]);

        if pi_idx == 1
            title(ax, caseNames{ci}, ...
                'Interpreter','latex', ...
                'FontSize',18, ...
                'FontWeight','normal');
        end

        if pci == 1
            ylabel(ax, samplePitchLabels{pi_idx}, ...
                'Interpreter','latex', ...
                'FontSize',17);
        else
            ylabel(ax, '');
        end

        xlabel(ax, '$x~(\mu\mathrm{m})$', ...
            'Interpreter','latex', ...
            'FontSize',13);

        set(ax, ...
            'FontSize',12, ...
            'TickLabelInterpreter','latex', ...
            'LineWidth',1.0, ...
            'TickDir','out');
end

end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 13;
cb.Label.String = 'Phase (rad)';
cb.Label.Interpreter = 'latex';

%% FIGURE 2 — PROPAGATED FOCAL ARRAYS

fig2 = figure('Color','w', 'Units','pixels', 'Position',[60 80 1350 1000]);

t2 = tiledlayout(n_pitch, n_plot_cases, ...
    'Padding','compact', ...
    'TileSpacing','compact');

xdraw = 90;

for pi_idx = 1:n_pitch

    for pci = 1:n_plot_cases

        ci = plotCaseIdx(pci);

        ax = nexttile;

        imagesc(ax, x_fiber*1e6, y_fiber*1e6, ...
            I_fiber_cases(:,:,pi_idx,ci));

        axis(ax, 'image');

        xlim(ax, [-xdraw xdraw]);
        ylim(ax, [-xdraw xdraw]);

        set(ax, 'YDir', 'normal');

        colormap(ax, turbo);
        caxis(ax, [0 1]);

        if pi_idx == 1
            title(ax, caseNames{ci}, ...
                'Interpreter','latex', ...
                'FontSize',18, ...
                'FontWeight','normal');
        end

        if pci == 1
            ylabel(ax, samplePitchLabels{pi_idx}, ...
                'Interpreter','latex', ...
                'FontSize',17);
        else
            ylabel(ax, '');
        end

        xlabel(ax, '$x~(\mu\mathrm{m})$', ...
            'Interpreter','latex', ...
            'FontSize',13);

        set(ax, ...
            'FontSize',12, ...
            'TickLabelInterpreter','latex', ...
            'LineWidth',1.0, ...
            'TickDir','out');
end

end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 13;
cb.Label.String = 'Normalized intensity';
cb.Label.Interpreter = 'latex';

%% FIGURE 3 — MODE OVERLAP EFFICIENCY

fig3 = figure('Color','w', 'Units','pixels', 'Position',[200 200 1200 800]);

t3 = tiledlayout(1,1, ...
    'Padding','compact', ...
    'TileSpacing','compact');

ax3plot = nexttile;

hold(ax3plot, 'on')

x_plot = [2 3 6 9 12 14];

for pi_idx = 1:n_pitch

    errorbar(ax3plot, x_plot, mean_overlap(pi_idx,:), std_overlap(pi_idx,:), ...
        'o-', ...
        'LineWidth',2, ...
        'MarkerSize',8, ...
        'CapSize',8);

end

hold(ax3plot, 'off')

grid(ax3plot, 'on')
box(ax3plot, 'on')

xlim(ax3plot, [1.5 14.8]);
ylim(ax3plot, [0 90]);

set(ax3plot, ...
    'XTick', x_plot, ...
    'XTickLabel', {'2','3','6','9','12','Continuous'}, ...
    'FontSize',20, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.2, ...
    'TickDir','out');

xlabel(ax3plot, 'Number of phase steps', ...
    'Interpreter','latex', ...
    'FontSize',22);

ylabel(ax3plot, 'Mean mode overlap efficiency (\%)', ...
    'Interpreter','latex', ...
    'FontSize',22);

legend(ax3plot, samplePitchLabels, ...
    'Interpreter','latex', ...
    'Location','northwest', ...
    'FontSize',20);

title(ax3plot, 'MS-MCF mode overlap after phase discretization and lateral sampling', ...
    'Interpreter','latex', ...
    'FontSize',24, ...
    'FontWeight','normal');
%% Representative phase masks and focal arrays
idx_pitch_3um   = find(abs(samplePitch_um_list - 3.0) < 1e-12, 1);
idx_pitch_400nm = find(abs(samplePitch_um_list - 0.4) < 1e-12, 1);

idx_case_3levels  = find(Nlevels_list == 3, 1);
idx_case_12levels = find(Nlevels_list == 12, 1);

fontLabel = 21;
fontTick  = 17;
fontCB    = 18;

xdraw_ms    = 90;
xdraw_fiber = 90;

rep_pitch_idx = [idx_pitch_3um, idx_pitch_3um, idx_pitch_400nm, idx_pitch_400nm];
rep_case_idx  = [idx_case_3levels, idx_case_12levels, idx_case_3levels, idx_case_12levels];

rep_names = { ...
    '3um_3levels', ...
    '3um_12levels', ...
    '400nm_3levels', ...
    '400nm_12levels'};

%% Phase masks

for ii = 1:numel(rep_names)

    pi_idx = rep_pitch_idx(ii);
    ci     = rep_case_idx(ii);

    fig_mask = figure('Color','w', 'Units','pixels', ...
        'Position',[100 100 620 560]);

    tmask = tiledlayout(1,1, ...
        'Padding','compact', ...
        'TileSpacing','compact');

    ax = nexttile;

    imagesc(ax, x_ms*1e6, y_ms*1e6, phase_masks(:,:,pi_idx,ci));

    axis(ax, 'image');
    xlim(ax, [-xdraw_ms xdraw_ms]);
    ylim(ax, [-xdraw_ms xdraw_ms]);
    set(ax, 'YDir', 'normal');

    colormap(ax, turbo);
    caxis(ax, [0 2*pi]);

    xlabel(ax, '$x~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',fontLabel);

    ylabel(ax, '$y~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',fontLabel);

    set(ax, ...
        'FontSize',fontTick, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1.1, ...
        'TickDir','out');
cb = colorbar(ax);
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fontCB;
    cb.Label.String = 'Phase (rad)';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = fontCB;

end

%% B — focal-array figures

for ii = 1:numel(rep_names)

    pi_idx = rep_pitch_idx(ii);
    ci     = rep_case_idx(ii);

    fig_focus = figure('Color','w', 'Units','pixels', ...
        'Position',[140 140 620 560]);

    tfocus = tiledlayout(1,1, ...
        'Padding','compact', ...
        'TileSpacing','compact');

    ax = nexttile;

    imagesc(ax, x_fiber*1e6, y_fiber*1e6, I_fiber_cases(:,:,pi_idx,ci));

    axis(ax, 'image');
    xlim(ax, [-xdraw_fiber xdraw_fiber]);
    ylim(ax, [-xdraw_fiber xdraw_fiber]);
    set(ax, 'YDir', 'normal');

    colormap(ax, turbo);
    caxis(ax, [0 1]);

    xlabel(ax, '$x~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',fontLabel);

    ylabel(ax, '$y~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',fontLabel);

    set(ax, ...
        'FontSize',fontTick, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1.1, ...
        'TickDir','out');
cb = colorbar(ax);
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = fontCB;
    cb.Label.String = 'Normalized intensity';
    cb.Label.Interpreter = 'latex';
    cb.Label.FontSize = fontCB;

end
