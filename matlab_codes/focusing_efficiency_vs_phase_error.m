% Effect of random inter-core phase errors on distal focusing efficiency.


clc
clear
close all

rng(1);

%% Fiber parameters

fiber_type = 'Spiral';

Lambda   = 6.7;       % [um]
CoreDiam = 2.5;       % [um]
lambda   = 0.920;     % [um]
n_cores  = 120;

n1 = 1.48;
n2 = 1.45;

taper = 0.67;
LambdaTaper = Lambda * taper;
r_core = CoreDiam * taper / 2;

%% Estimate mode size

k0 = 2*pi/lambda;
V = k0 * r_core * sqrt(n1^2 - n2^2);
w0 = r_core * (0.65 + 1.619/V^(3/2) + 2.879/V^6);

%% Generate MCF geometry

fiber = FiberClassGenerator(fiber_type, n_cores, LambdaTaper);

%% Simulation grid

z_focus = 400;       % [um]

N = 1000;

x = linspace(-200, 200, N);
y = x;

dx = abs(x(2) - x(1));
dy = dx;

xf = linspace(-1/(2*dx), 1/(2*dx), N);
yf = xf;

[X, Y] = meshgrid(x, y);
R = sqrt(X.^2 + Y.^2);

k = 2*pi/lambda;

[Xf, Yf] = meshgrid(lambda*z_focus*xf, lambda*z_focus*yf);

%% Ideal focusing phase profile

target_pos = [0 0 z_focus];

phi_focus = zeros(n_cores, 1);

for i = 1:n_cores

    phi_focus(i) = CalculatePhaseCore( ...
        fiber.posarray(i,1), ...
        fiber.posarray(i,2), ...
        target_pos, ...
        lambda);

end

%% Precompute Gaussian core modes

CoreModes = zeros(N*N, n_cores, 'single');

for i = 1:n_cores

    R_core = sqrt( ...
        (X - fiber.posarray(i,1)).^2 + ...
        (Y - fiber.posarray(i,2)).^2 );

    mode_i = fiber.A(i) .* exp(-(R_core).^2 / w0^2);

    CoreModes(:,i) = single(mode_i(:));

end

%% Phase-error levels

sigma_phi_list = 0:0.05:1.0;       % [rad]
n_cases = numel(sigma_phi_list);

sigma_examples = [0 0.25 0.50 1.00];   % [rad]

exampleNames = { ...
    'Ideal', ...
    '$\sigma_\phi = 0.25~\mathrm{rad}$', ...
    '$\sigma_\phi = 0.50~\mathrm{rad}$', ...
    '$\sigma_\phi = 1.00~\mathrm{rad}$'};

n_examples = numel(sigma_examples);

nReal = 10;

Strehl_all = zeros(nReal, n_cases);
relative_Strehl_all = zeros(nReal, n_cases);

IntensityExample = zeros(N, N, n_examples, 'single');

%% Simulation loop

for rr = 1:nReal

    % Random piston phase profile over the 120 cores
    random_piston = randn(n_cores, 1);

    % Remove global piston and normalize RMS to 1 rad
    random_piston = random_piston - mean(random_piston);
    random_piston = random_piston ./ sqrt(mean(random_piston.^2));

    for c = 1:n_cases

        sigma_phi = sigma_phi_list(c);

        % Injected random phase error in rad
        phase_error = sigma_phi * random_piston;

        % Focusing phase + fabrication phase error
        phi_test = phi_focus + phase_error;

        weights = exp(1i * phi_test);

        E_vec = complex(zeros(N*N,1,'single'));

        for i = 1:n_cores

            E_vec = E_vec + CoreModes(:,i) .* single(weights(i));

        end

        E_fiber = reshape(E_vec, N, N);

        % Fresnel propagation
        E_fresnel = E_fiber .* exp(1i * R.^2 * k / (2*z_focus));

        E_propag = ...
            dx * dy * ...
            fftshift(fft2(ifftshift(E_fresnel))) * ...
            exp(1i*k*z_focus) / (1i*lambda*z_focus) .* ...
            exp(1i*(Xf.^2 + Yf.^2) * k / (2*z_focus));

        I = abs(E_propag).^2;

        [Strehl_tmp, ~, ~, ~] = GetStrehl_useMax_Generic(double(I), 0);

        Strehl_all(rr,c) = Strehl_tmp * 100;

        % Store focal spots only for the first repetition
        if rr == 1

            idxExample = find(abs(sigma_phi - sigma_examples) < 1e-12);

            if ~isempty(idxExample)

                IntensityExample(:,:,idxExample) = single(I);

            end

        end

    end

    relative_Strehl_all(rr,:) = Strehl_all(rr,:) ./ Strehl_all(rr,1);

end

mean_relative_Strehl = mean(relative_Strehl_all, 1);
std_relative_Strehl  = std(relative_Strehl_all, [], 1);

mean_Strehl = mean(Strehl_all, 1);
std_Strehl  = std(Strehl_all, [], 1);

%% Figure 1 — focal spots

fig1 = figure('Color','w', 'Units','pixels', 'Position',[80 120 1250 360]);

tiledlayout(1, n_examples, ...
    'Padding','compact', ...
    'TileSpacing','compact');

xdraw = 60;

I_ref_max = max(IntensityExample(:,:,1), [], 'all');

for e = 1:n_examples

    ax = nexttile;

    Iplot = double(IntensityExample(:,:,e)) ./ I_ref_max;

    imagesc(ax, Xf(1,:), Yf(:,1), Iplot);

    axis(ax, 'image');

    xlim(ax, [-xdraw xdraw]);
    ylim(ax, [-xdraw xdraw]);

    set(ax, 'YDir', 'normal');

    colormap(ax, turbo);
    caxis(ax, [0 1]);

    title(ax, exampleNames{e}, ...
        'Interpreter','latex', ...
        'FontSize',18, ...
        'FontWeight','normal');

    xlabel(ax, '$x~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',16);

    ylabel(ax, '$y~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',16);

    set(ax, ...
        'FontSize',14, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1.0, ...
        'TickDir','out');

end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 14;
cb.Label.String = 'Normalized intensity';
cb.Label.Interpreter = 'latex';

%% Figure 2 — relative focusing efficiency

fig2 = figure('Color','w', 'Units','pixels', 'Position',[180 170 850 520]);

errorbar(sigma_phi_list, mean_relative_Strehl, std_relative_Strehl, ...
    'o-', ...
    'LineWidth',2, ...
    'MarkerSize',7, ...
    'CapSize',6);

grid on
box on

xlim([0 max(sigma_phi_list)]);
ylim([0 1.05]);

xlabel('Injected phase error $\sigma_\phi$ (rad)', ...
    'Interpreter','latex', ...
    'FontSize',20);

ylabel('Relative focusing efficiency', ...
    'Interpreter','latex', ...
    'FontSize',20);

title('Relative focusing efficiency loss from random phase errors', ...
    'Interpreter','latex', ...
    'FontSize',20, ...
    'FontWeight','normal');

set(gca, ...
    'FontSize',18, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.2, ...
    'TickDir','out');
