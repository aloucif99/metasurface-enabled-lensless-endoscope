% Mode-overlap efficiency versus axial propagation distance.
% local_voronoi_segmentation.m in the same folder.

clear
close all
clc

%% Parameters

lambda = 920e-9;          % [m]
period_MS = 400e-9;       % [m]
n_cores = 120;

size_plane = 180e-6;      % [m]
N = round(size_plane / period_MS);

x = period_MS * ((-N/2):(N/2-1));
[X, Y] = meshgrid(x, x);

z0_MS = -30e-6;           % [m]
phi2_MS = pi/(lambda*z0_MS) * period_MS^2;

z = 10e-6:2e-6:70e-6;    % [m]

%% Calibrated inter-core phases

load(fullfile(fileparts(mfilename('fullpath')), 'PhaseVector400um.mat'),'phase_vec');
Phi = phase_vec(:);
Phi = Phi(1:n_cores);

%% Fermat core positions

p_MS = 30;
size_beamlet_ms = p_MS / 1.72;

rho = sqrt(1:1000)';
theta = (3-sqrt(5)) * pi * (1:1000)';
poscore = size_beamlet_ms * [rho.*cos(theta), rho.*sin(theta)];

%% Metasurface phase mask

basevecX = exp(1i*Phi);

phase_mask = phase_mask_mcf( basevecX, poscore, phi2_MS, N, N, 0);

phase_mask = squeeze(phase_mask);

%% Fiber mode - Parabolic GRIN fiber

a = 1.25e-6;              % Core radius [m]
n_core = 1.48;
n_clad = 1.45;

g_grin = 2;               % Parabolic GRIN profile

V = 2*pi*a/lambda * sqrt(n_core^2 - n_clad^2);

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
%% Voronoi regions around the MCF cores

[IN, SegSizes] = local_voronoi_segmentation(poscore, N, N);
Npix = N^2;

%% Propagation

fx = 1/(N*period_MS) * ((-N/2):(N/2-1));
[Fx, Fy] = meshgrid(fx, fx);

field_fft = fftshift(fft2(ifftshift(exp(1i*phase_mask))));

overlap = zeros(n_cores, numel(z));

for iz = 1:numel(z)

    H = exp(-1i*pi*lambda*z(iz) * (Fx.^2 + Fy.^2));

    field = fftshift(ifft2(ifftshift(field_fft .* H)));
    field = field / sqrt(sum(abs(field(:)).^2));

    field_fft_fiber = fftshift(fft2(ifftshift(field)));
    C = fftshift(ifft2(ifftshift(fiber_mode_fft .* field_fft_fiber)));

    for k = 1:n_cores
        idx = IN{k};
        overlap(k,iz) = max(abs(C(idx)).^2) * ...
            (Npix / SegSizes(k));
    end
end

overlap = 100 * overlap;

%% Mean overlap versus axial distance

mean_overlap = mean(overlap,1);
std_overlap = std(overlap,[],1);

[~,idx_best] = max(mean_overlap);

figure
errorbar(z*1e6, mean_overlap, std_overlap, '-o', ...
    'LineWidth',1.5,'MarkerSize',6)
xlabel('Observation distance z (\mum)')
ylabel('Mode-overlap efficiency (%)')
ylim([0 100])
grid on

%% Distribution across the 120 cores at the best plane

figure
histogram(overlap(:,idx_best),'BinWidth',2)
xlabel('Mode-overlap efficiency (%)')
ylabel('Number of cores')
grid on
