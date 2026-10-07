% Lateral MS-MCF misalignment: core maps and mean mode-overlap efficiency.
% Requires PhaseVector400um.mat, phase_mask_mcf.m and
% local_voronoi_segmentation.m in the same folder.

clear
close all
clc


%% Parameters

lambda = 920e-9;
period_MS = 400e-9;
n_cores = 120;

size_plane = 180e-6;
N = round(size_plane/period_MS);

x = period_MS*((-N/2):(N/2-1));
[X,Y] = meshgrid(x,x);

z0_MS = -30e-6;
f_fixed = 28e-6;

x_offset_um = 0:0.25:3;
x_offset_list = x_offset_um*1e-6;

display_offsets_um = [0 0.5 1 2];

%% Phase vector and Fermat geometry

phasevec_path = fullfile(fileparts(mfilename('fullpath')), 'PhaseVector400um.mat');

load(phasevec_path,'phase_vec');
Phi = phase_vec(:);
Phi = Phi(1:n_cores);

p_MS = 30;
size_beamlet = p_MS/1.72;

rho = sqrt(1:1000)';
theta = (3-sqrt(5))*pi*(1:1000)';
poscore = size_beamlet*[rho.*cos(theta), rho.*sin(theta)];

basevecX = exp(1i*Phi);

%% Fiber mode - Parabolic GRIN fiber

a = 1.25e-6;              % Core radius [m]
n_core = 1.48;
n_clad = 1.45;

g_grin = 2;               % Parabolic GRIN profile
Delta_grin = (n_core^2 - n_clad^2) / (2*n_core^2);
V = 2*pi/lambda * n_core * a * sqrt(2*Delta_grin);

% Marcuse Gaussian approximation for graded-index fiber
A_grin = sqrt(2/5 * (1 + 4*(2/g_grin)^(5/6)));
B_grin = exp(0.298/g_grin) - 1 ...
       + 1.478*(1 - exp(-0.077*g_grin));
C_grin = 3.76 + exp(4.19/g_grin^0.418);

w_fiber = a * ( ...
    A_grin / V^(2/(g_grin+2)) + ...
    B_grin / V^(3/2) + ...
    C_grin / V^6 );

fiber_mode = exp(-(X.^2 + Y.^2)/w_fiber^2);
fiber_mode = fiber_mode / sqrt(sum(abs(fiber_mode(:)).^2));

fiber_mode_fft = fftshift(fft2(ifftshift(fiber_mode)));

%% Core regions

poscore_phys = poscore*period_MS;
poscore_pix = poscore_phys/period_MS;

[~,SegSizes] = local_voronoi_segmentation(poscore_pix,N,N);
Npix = N^2;

%% Metasurface field at the coupling plane

phi2 = pi/(lambda*z0_MS)*period_MS^2;

phase_mask = phase_mask_mcf(basevecX,poscore,phi2,N,N,0);
phase_mask = squeeze(phase_mask);

fx = 1/(N*period_MS)*((-N/2):(N/2-1));
[Fx,Fy] = meshgrid(fx,fx);

field_fft = fftshift(fft2(ifftshift(exp(1i*phase_mask))));
H = exp(-1i*pi*lambda*f_fixed*(Fx.^2+Fy.^2));

spotfield = fftshift(ifft2(ifftshift(field_fft.*H)));

%% Lateral-offset sweep

results_overlap = NaN(n_cores,numel(x_offset_list));
mean_overlap = NaN(size(x_offset_list));
std_overlap = NaN(size(x_offset_list));

coupling_maps = NaN(n_cores,numel(display_offsets_um));

for i = 1:numel(x_offset_list)

    x_offset = x_offset_list(i);

    shifted_field = interp2( ...
        X,Y,spotfield,X-x_offset,Y,'linear',0);

    shifted_field = shifted_field / ...
        sqrt(sum(abs(shifted_field(:)).^2));

    S_fft = fftshift(fft2(ifftshift(shifted_field)));
    overlap_map = fftshift(ifft2(ifftshift( ...
        fiber_mode_fft.*S_fft)));

    coupleff = NaN(n_cores,1);

    for k = 1:n_cores

        x_core = poscore_phys(k,1);
        y_core = poscore_phys(k,2);

        [~,col] = min(abs(x-x_core));
        [~,row] = min(abs(x-y_core));

        coupleff(k) = abs(overlap_map(row,col)).^2 * ...
            (Npix/SegSizes(k));
    end

    results_overlap(:,i) = 100*coupleff;
    mean_overlap(i) = mean(results_overlap(:,i),'omitnan');
    std_overlap(i) = std(results_overlap(:,i),[],'omitnan');

    j = find(abs(display_offsets_um-x_offset_um(i))<1e-12);
    if ~isempty(j)
        coupling_maps(:,j) = results_overlap(:,i);
    end
end

%% Core maps at selected offsets

x_core_um = poscore_phys(1:n_cores,1)*1e6;
y_core_um = poscore_phys(1:n_cores,2)*1e6;

fig1 = figure('Color','w','Position',[80 80 800 800]);

tiledlayout(2,2, ...
    'Padding','compact', ...
    'TileSpacing','compact');

for j = 1:numel(display_offsets_um)

    ax = nexttile;

    scatter(ax, x_core_um, y_core_um, 95, coupling_maps(:,j), 'filled');

    axis(ax,'image');
    xlim(ax,[-90 90]);
    ylim(ax,[-90 90]);

    colormap(ax,"parula");
    caxis(ax,[0 100]);

    xlabel(ax,'$x~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',20);

    ylabel(ax,'$y~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',20);

    title(ax,sprintf('$x_{\\mathrm{offset}} = %g~\\mu\\mathrm{m}$', ...
        display_offsets_um(j)), ...
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

%% Mean overlap versus lateral offset

fig2 = figure('Color','w','Position',[120 120 980 530]);

errorbar(x_offset_um, mean_overlap, std_overlap, ...
    '-o', ...
    'LineWidth',2, ...
    'MarkerSize',7, ...
    'CapSize',7);

grid on
box on

xlabel('Lateral offset $x_{\mathrm{offset}}~(\mu\mathrm{m})$', ...
    'Interpreter','latex', ...
    'FontSize',20);

ylabel('Mean overlap efficiency (\%)', ...
    'Interpreter','latex', ...
    'FontSize',20);

title('Overlap efficiency vs lateral MS - MCF misalignment', ...
    'Interpreter','latex', ...
    'FontSize',23, ...
    'FontWeight','normal');

legend('Mean $\pm$ std over 120 cores', ...
    'Interpreter','latex', ...
    'Location','best', ...
    'FontSize',20);

xlim([0 3]);
xticks(0:0.5:3);
ylim([0 100]);

set(gca, ...
    'FontSize',22, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.2, ...
    'TickDir','out');
