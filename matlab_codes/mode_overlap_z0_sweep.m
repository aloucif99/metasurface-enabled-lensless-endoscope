% Mode-overlap efficiency versus MS design focal distance.


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

z0_list = -20e-6:-5e-6:-50e-6;   % [m]
z = 10e-6:2e-6:70e-6;             % [m]

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

%% Voronoi regions around the MCF cores

[IN, SegSizes] = local_voronoi_segmentation(poscore, N, N);
Npix = N^2;

%% Fourier grid

fx = 1/(N*period_MS) * ((-N/2):(N/2-1));
[Fx, Fy] = meshgrid(fx, fx);

%% Sweep over MS design focal distance

overlap_all = zeros(n_cores, numel(z), numel(z0_list));
mean_overlap = zeros(numel(z0_list), numel(z));
std_overlap = zeros(numel(z0_list), numel(z));

for iz0 = 1:numel(z0_list)

    z0 = z0_list(iz0);
    phi2 = pi/(lambda*z0) * period_MS^2;

    phase_mask = phase_mask_mcf( ...
        basevecX, poscore, phi2, N, N, 0);

    phase_mask = squeeze(phase_mask);

    field_fft = fftshift(fft2(ifftshift(exp(1i*phase_mask))));

    for iz = 1:numel(z)

        H = exp(-1i*pi*lambda*z(iz) * (Fx.^2 + Fy.^2));

        field = fftshift(ifft2(ifftshift(field_fft .* H)));
        field = field / sqrt(sum(abs(field(:)).^2));

        field_fft_fiber = fftshift(fft2(ifftshift(field)));
        C = fftshift(ifft2(ifftshift( ...
            fiber_mode_fft .* field_fft_fiber)));

        for k = 1:n_cores
            idx = IN{k};
            overlap_all(k,iz,iz0) = ...
                max(abs(C(idx)).^2) * (Npix / SegSizes(k));
        end
    end

    overlap_all(:,:,iz0) = 100 * overlap_all(:,:,iz0);
    mean_overlap(iz0,:) = mean(overlap_all(:,:,iz0),1);
    std_overlap(iz0,:) = std(overlap_all(:,:,iz0),[],1);
end

%% Best overlap for each design focal distance

best_mean = zeros(size(z0_list));
best_std = zeros(size(z0_list));
best_z = zeros(size(z0_list));
best_idx = zeros(size(z0_list));

for iz0 = 1:numel(z0_list)
    [best_mean(iz0),best_idx(iz0)] = max(mean_overlap(iz0,:));
    best_std(iz0) = std_overlap(iz0,best_idx(iz0));
    best_z(iz0) = z(best_idx(iz0));
end

figure
errorbar(z0_list*1e6, best_mean, best_std, '-o', ...
    'LineWidth',1.5,'MarkerSize',6)
xlabel('Design focal distance z_0^{MS} (\mum)')
ylabel('Best mean mode-overlap efficiency (%)')
ylim([0 100])
grid on

%% Axial overlap curve for z0 = -30 um

[~,idx_30] = min(abs(z0_list + 30e-6));

figure
errorbar(z*1e6, mean_overlap(idx_30,:), std_overlap(idx_30,:), ...
    '-o','LineWidth',1.5,'MarkerSize',5)
xlabel('Observation distance z (\mum)')
ylabel('Mode-overlap efficiency (%)')
ylim([0 100])
grid on

%% Distribution across the 120 cores at the best plane for z0 = -30 um

data_best = overlap_all(:,best_idx(idx_30),idx_30);

figure
histogram(data_best,'BinWidth',2)
xlabel('Mode-overlap efficiency (%)')
ylabel('Number of cores')
grid on
