% Interpolate the simulated GaN nanopillar transmission and phase lookup table.
% Requires pillars_GaNonSaph_wv920_prd400_thk1500nm.mat in the same folder.

clear; close all; clc;

%% Load data

load(fullfile(fileparts(mfilename('fullpath')), 'pillars_GaNonSaph_wv920_prd400_thk1500nm.mat'));

% Some versions of the simulation file use DeltaXi for the phase array.
if exist('DeltaXi','var') && ~exist('DeltaX','var')
    DeltaX = DeltaXi;
end

R = R(:).';
z_max_p = z_max_p(:).';

if size(Tr,1) ~= numel(R)
    Tr = Tr.';
end
if size(DeltaX,1) ~= numel(R)
    DeltaX = DeltaX.';
end

%% Phase processing

DeltaX_unwrapped = unwrap(DeltaX, [], 1);

ref_phase  = DeltaX_unwrapped(1,1);
DeltaX_ref = DeltaX_unwrapped - ref_phase;

%% Interpolation

[R_grid, H_grid] = ndgrid(R, z_max_p);

R_fine = linspace(min(R), max(R), 2000);
H_fine = linspace(min(z_max_p), max(z_max_p), 2000);

[R_fine_grid, H_fine_grid] = ndgrid(R_fine, H_fine);

F_Tr = griddedInterpolant(R_grid, H_grid, Tr, 'spline');
F_Phase = griddedInterpolant(R_grid, H_grid, DeltaX_ref, 'spline');

Tr_fine = F_Tr(R_fine_grid, H_fine_grid);
Phase_fine_unwrapped = F_Phase(R_fine_grid, H_fine_grid);
Phase_fine_wrapped = mod(Phase_fine_unwrapped, 2*pi);

%% Display range and selected values

D_fine_nm = 2 * R_fine * 1e9;
H_fine_nm = H_fine * 1e9;

Dmin_nm = 100;
Dmax_nm = 300;

idxD = (D_fine_nm >= Dmin_nm) & (D_fine_nm <= Dmax_nm);

D_sel_nm = linspace(100, 300, 12);
R_sel_m  = (D_sel_nm * 1e-9) / 2;

H_target_nm = 1500;
H_target_m  = H_target_nm * 1e-9;

Tr_cut = F_Tr(R_fine, H_target_m * ones(size(R_fine)));
Phase_cut_wrap = mod(F_Phase(R_fine, H_target_m * ones(size(R_fine))), 2*pi);

Tr_sel = F_Tr(R_sel_m, H_target_m * ones(size(R_sel_m)));
Phase_sel_wrap = mod(F_Phase(R_sel_m, H_target_m * ones(size(R_sel_m))), 2*pi);

%% Plot style

fontTitle = 24;
fontLabel = 24;
fontTick  = 21;
fontCB    = 21;

lw_curve = 2.4;
lw_marker = 2.0;
ms_sel = 8;

figW_map = 18;
figH_map = 14;

figW_cut = 18;
figH_cut = 14;

%% Transmission map

fig1 = figure('Color','w', 'Units','centimeters', ...
    'Position',[2 2 figW_map figH_map]);

t1 = tiledlayout(fig1, 1, 1, ...
    'Padding','compact', ...
    'TileSpacing','compact');

ax1 = nexttile;

imagesc(ax1, D_fine_nm, H_fine_nm/1000, Tr_fine.');
axis(ax1, 'xy');

xlim(ax1, [Dmin_nm Dmax_nm]);
ylim(ax1, [min(H_fine_nm) max(H_fine_nm)]/1000);
caxis(ax1, [0 1]);

colormap(ax1, parula);

cb1 = colorbar(ax1);
cb1.TickLabelInterpreter = 'latex';
cb1.FontSize = fontCB;
cb1.Label.String = 'Transmission';
cb1.Label.Interpreter = 'latex';
cb1.Label.FontSize = fontCB;

xlabel(ax1, 'Diameter (nm)', ...
    'Interpreter','latex', ...
    'FontSize',fontLabel);

ylabel(ax1, 'Height ($\mu$m)', ...
    'Interpreter','latex', ...
    'FontSize',fontLabel);

title(ax1, 'Transmission', ...
    'Interpreter','latex', ...
    'FontSize',fontTitle, ...
    'FontWeight','normal');

set(ax1, ...
    'FontSize',fontTick, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.3, ...
    'TickDir','out', ...
    'Box','on');
%% Phase map

fig2 = figure('Color','w', 'Units','centimeters', ...
    'Position',[2 2 figW_map figH_map]);

t2 = tiledlayout(fig2, 1, 1, ...
    'Padding','compact', ...
    'TileSpacing','compact');

ax2 = nexttile;

imagesc(ax2, D_fine_nm, H_fine_nm/1000, Phase_fine_wrapped.');
axis(ax2, 'xy');

xlim(ax2, [Dmin_nm Dmax_nm]);
ylim(ax2, [min(H_fine_nm) max(H_fine_nm)]/1000);
caxis(ax2, [0 2*pi]);

colormap(ax2, turbo);

cb2 = colorbar(ax2);
cb2.TickLabelInterpreter = 'latex';
cb2.FontSize = fontCB;
cb2.Ticks = [0 pi/2 pi 3*pi/2 2*pi];
cb2.TickLabels = {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'};
cb2.Label.String = 'Phase (rad)';
cb2.Label.Interpreter = 'latex';
cb2.Label.FontSize = fontCB;

xlabel(ax2, 'Diameter (nm)', ...
    'Interpreter','latex', ...
    'FontSize',fontLabel);

ylabel(ax2, 'Height ($\mu$m)', ...
    'Interpreter','latex', ...
    'FontSize',fontLabel);

title(ax2, 'Phase', ...
    'Interpreter','latex', ...
    'FontSize',fontTitle, ...
    'FontWeight','normal');

set(ax2, ...
    'FontSize',fontTick, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.3, ...
    'TickDir','out', ...
    'Box','on');
%% Transmission at H = 1.5 um

fig3 = figure('Color','w', 'Units','centimeters', ...
    'Position',[2 2 figW_cut figH_cut]);

t3 = tiledlayout(fig3, 1, 1, ...
    'Padding','compact', ...
    'TileSpacing','compact');

ax3 = nexttile;

plot(ax3, D_fine_nm(idxD), Tr_cut(idxD), '-', ...
    'LineWidth', lw_curve);

hold(ax3, 'on');

plot(ax3, D_sel_nm, Tr_sel, 'o', ...
    'LineWidth', lw_marker, ...
    'MarkerSize', ms_sel, ...
    'MarkerFaceColor', 'none');

hold(ax3, 'off');

grid(ax3, 'on');
box(ax3, 'on');

xlim(ax3, [Dmin_nm Dmax_nm]);
ylim(ax3, [0 1.05]);

xlabel(ax3, 'Diameter (nm)', ...
    'Interpreter','latex', ...
    'FontSize',fontLabel);

ylabel(ax3, 'Transmission', ...
    'Interpreter','latex', ...
    'FontSize',fontLabel);

legend(ax3, {'Simulation', 'Selected diameters'}, ...
    'Interpreter','latex', ...
    'Location','southwest', ...
    'FontSize',18, ...
    'Box','off');

set(ax3, ...
    'FontSize',fontTick, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.3, ...
    'TickDir','out');
%% Phase at H = 1.5 um

fig4 = figure('Color','w', 'Units','centimeters', ...
    'Position',[2 2 figW_cut figH_cut]);

t4 = tiledlayout(fig4, 1, 1, ...
    'Padding','compact', ...
    'TileSpacing','compact');

ax4 = nexttile;

plot(ax4, D_fine_nm(idxD), Phase_cut_wrap(idxD), '-', ...
    'LineWidth', lw_curve);

hold(ax4, 'on');

plot(ax4, D_sel_nm, Phase_sel_wrap, 'o', ...
    'LineWidth', lw_marker, ...
    'MarkerSize', ms_sel, ...
    'MarkerFaceColor', 'none');

hold(ax4, 'off');

grid(ax4, 'on');
box(ax4, 'on');

xlim(ax4, [Dmin_nm Dmax_nm]);
ylim(ax4, [0 2*pi]);

yticks(ax4, [0 pi/2 pi 3*pi/2 2*pi]);
yticklabels(ax4, {'$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'});

xlabel(ax4, 'Diameter (nm)', ...
    'Interpreter','latex', ...
    'FontSize',fontLabel);

ylabel(ax4, 'Phase (rad)', ...
    'Interpreter','latex', ...
    'FontSize',fontLabel);

legend(ax4, {'Simulation', 'Selected diameters'}, ...
    'Interpreter','latex', ...
    'Location','northwest', ...
    'FontSize',18, ...
    'Box','off');

set(ax4, ...
    'FontSize',fontTick, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.3, ...
    'TickDir','out');
