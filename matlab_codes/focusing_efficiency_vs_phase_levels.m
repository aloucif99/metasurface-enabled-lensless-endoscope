% Effect of phase discretization on distal focusing efficiency.
% The random-number generator is fixed for reproducibility.

clc
clear
close all

rng(1);

%% Fiber parameters

fiber_type = 'Spiral';

Lambda   = 6.7;      % [um]
CoreDiam = 2.5;      % [um]
lambda   = 0.920;    % [um]
n_cores  = 120;

n1 = 1.48;
n2 = 1.45;

taper = 0.67;
LambdaTaper = Lambda * taper;
r_core = CoreDiam * taper / 2;

%% Mode size estimate

k0 = 2*pi/lambda;
V  = k0 * r_core * sqrt(n1^2 - n2^2);

w0 = r_core * ...
    (0.65 + 1.619/(V^(3/2)) + 2.879/(V^6));

%% MCF geometry

fiber = FiberClassGenerator(fiber_type, n_cores, LambdaTaper);

%% Simulation grid

z_focus = 400;      % [um]

N = 1000;

x = linspace(-200, 200, N);
y = x;

dx = abs(x(2)-x(1));
dy = dx;

xf = linspace(-1/(2*dx), 1/(2*dx), N);
yf = xf;

[X, Y] = meshgrid(x, y);
R = sqrt(X.^2 + Y.^2);

k = 2*pi/lambda;

[Xf, Yf] = meshgrid(lambda*z_focus*xf, lambda*z_focus*yf);

%% Ideal focusing phase

target_pos = [0 0 z_focus];

phi_focus = zeros(n_cores,1);

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

%% Test cases

caseNames = {'2 levels', '3 levels', '6 levels', '9 levels', '12 levels', 'Continuous'};

Nlevels_list = [2 3 6 9 12 Inf];

n_cases = numel(Nlevels_list);
nReal = 10;

Strehl_all = zeros(nReal, n_cases);
relative_Strehl_all = zeros(nReal, n_cases);

IntensityExample = zeros(N, N, n_cases, 'single');

%% Main loop

for rr = 1:nReal

    % Random phase uniformly distributed in [0, 2pi)
    phi_random = 2*pi*rand(n_cores,1);

    for c = 1:n_cases

        if isinf(Nlevels_list(c))

            % Ideal continuous focusing phase, no random phase
            phi_actual = phi_focus;

        else

            Nlevels = Nlevels_list(c);
            dphi = 2*pi / Nlevels;

            % Required phase compensation:
            % after random phase + compensation, we want phi_focus
            phi_required = phi_focus - phi_random;

            % Wrap to [0, 2pi)
            phi_required_wrapped = mod(phi_required, 2*pi);

            % Quantize to nearest available phase level
            level_index = round(phi_required_wrapped / dphi);

            % 2pi equivalent to 0
            level_index = mod(level_index, Nlevels);

            phi_comp = level_index * dphi;

            % Actual phase after random phase + discretized compensation
            phi_actual = phi_random + phi_comp;

        end

        %% -----------------------------------------------------------------
        % Build distal field
        % -----------------------------------------------------------------

        weights = exp(1i * phi_actual);

        E_vec = complex(zeros(N*N,1,'single'));

        for i = 1:n_cores

            E_vec = E_vec + CoreModes(:,i) .* single(weights(i));

        end

        E_fiber = reshape(E_vec, N, N);

        %% -----------------------------------------------------------------
        % Fresnel propagation
        % -----------------------------------------------------------------

        E_fresnel = E_fiber .* exp(1i * R.^2 * k / (2*z_focus));

        E_propag = ...
            dx * dy * ...
            fftshift(fft2(ifftshift(E_fresnel))) * ...
            exp(1i*k*z_focus) / (1i*lambda*z_focus) .* ...
            exp(1i*(Xf.^2 + Yf.^2) * k / (2*z_focus));

        I = abs(E_propag).^2;

        [Strehl_tmp, ~, ~, ~] = GetStrehl_useMax_Generic(double(I), 0);

        Strehl_all(rr,c) = Strehl_tmp * 100;

        % Store focal spots from first repetition only
        if rr == 1

            IntensityExample(:,:,c) = single(I);

        end

    end

    % Normalize all cases by the ideal case, which is the last one
    relative_Strehl_all(rr,:) = Strehl_all(rr,:) ./ Strehl_all(rr,end);

end

%% Statistics

mean_Strehl = mean(Strehl_all, 1);
std_Strehl  = std(Strehl_all, 0, 1);

mean_relative_Strehl = mean(relative_Strehl_all, 1);
std_relative_Strehl  = std(relative_Strehl_all, 0, 1);

%% Figure 1 — focal spots

fig1 = figure('Color','w', 'Units','pixels', 'Position',[2 2 1000 700]);

t1 = tiledlayout(2, 3, ...
    'Padding','compact', ...
    'TileSpacing','compact');

xdraw = 60;

I_ref_max = max(IntensityExample(:,:,end), [], 'all');

fontTitle = 20;
fontLabel = 19;
fontTick  = 17;
fontCB    = 20;

for c = 1:n_cases

    ax = nexttile;

    Iplot = double(IntensityExample(:,:,c)) ./ I_ref_max;

    imagesc(ax, Xf(1,:), Yf(:,1), Iplot);

    axis(ax, 'image');

    xlim(ax, [-xdraw xdraw]);
    ylim(ax, [-xdraw xdraw]);

    set(ax, 'YDir', 'normal');

    colormap(ax, turbo);
    caxis(ax, [0 1]);

    title(ax, caseNames{c}, ...
        'Interpreter','latex', ...
        'FontSize',fontTitle, ...
        'FontWeight','normal');

    xlabel(ax, '$x~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',fontLabel);

    ylabel(ax, '$y~(\mu\mathrm{m})$', ...
        'Interpreter','latex', ...
        'FontSize',fontLabel);

    set(ax, ...
        'FontSize',fontTick, ...
        'TickLabelInterpreter','latex', ...
        'LineWidth',1.0, ...
        'TickDir','out');
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = fontCB;
cb.Label.String = 'Normalized intensity';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = fontCB;

%% Figure 2 — relative focusing efficiency

fig2 = figure('Color','w', 'Units','pixels', 'Position',[50 50 950  450]);

t2 = tiledlayout(1,1, ...
    'Padding','compact', ...
    'TileSpacing','compact');

ax2 = nexttile;

% Real x positions for discrete cases + one attached point for ideal
x_plot = [2 3 6 9 12 14];

errorbar(ax2, x_plot, mean_relative_Strehl, std_relative_Strehl, ...
    'o-', ...
    'LineWidth',2, ...
    'MarkerSize',8, ...
    'CapSize',8);

grid(ax2, 'on')
box(ax2, 'on')

xlim(ax2, [1.5 14.8]);
ylim(ax2, [0 1.1]);

set(ax2, ...
    'XTick', x_plot, ...
    'XTickLabel', {'2','3','6','9','12','Ideal'}, ...
    'FontSize',18, ...
    'TickLabelInterpreter','latex', ...
    'LineWidth',1.2, ...
    'TickDir','out');

xlabel(ax2, 'Number of phase steps', ...
    'Interpreter','latex', ...
    'FontSize',22);

ylabel(ax2, 'Relative focusing efficiency', ...
    'Interpreter','latex', ...
    'FontSize',22);

title(ax2, 'Effect of discrete phase compensation on focusing efficiency', ...
    'Interpreter','latex', ...
    'FontSize',22, ...
    'FontWeight','normal');
