clc
clear all
close all
%addpath('C:\Luca\MatlabCodes\PropagationSimulation\FunctionsSimulation');
z_focus = 400;
%%%%%%%%%%%%%%%%%%%
%%% FIBER PARAM %%%
%%%%%%%%%%%%%%%%%%%
% % Fiber used for Genchi & Hofer, 2025
% fiber_type = 'Spiral';
% Lambda = 10.5; % um
% CoreDiam = 5; % um
% lambda = 0.920; % um
% n_cores = 120;
% n1 = 1.48;
% n2 = 1.45;
% taper = 0.40;

% Fiber used for Loucif & Genchi, 2026
fiber_type = 'Spiral';
Lambda = 6.7;
CoreDiam = 2.5;
lambda = 0.920;
n_cores = 120;
n1 = 1.48;
n2 = 1.45;
% taper = 1;
taper = 0.66;

StepIndex = 0; % 0 = GRIN, 1 = Step Index
deltan = n1-n2;
CoreNA = sqrt(n1^2 - n2^2);
msg_setup = sprintf(strcat(fiber_type, sprintf(' fiber parameters:\nLambda = %.1f um, core diam. = %.1f um, n1 = %.2f, n2 = %.2f, theoretical NA = %.2f.', Lambda, CoreDiam, n1, n2, CoreNA)));

%%%%%%%%%%%%%
%%% TAPER %%%
%%%%%%%%%%%%%
if taper >= 1
    LambdaTaper = Lambda;
    CoreDiamOutput = CoreDiam;
    r_core = CoreDiam / 2;
    TaperString = 'NoTaper';
    msg_taper = sprintf('No taper for this fiber.')
else
    LambdaTaper = Lambda * taper;
    CoreDiamOutput = CoreDiam * taper;
    r_core = CoreDiam * taper / 2; % radius of the core
    TaperString = strcat('Taper0',num2str(taper*100,'%d'));
    msg_taper = sprintf('Taper of %.2f selected. Therefore at the output:\nLambda = %.1f um, core diam. = %.2f um\n\n', taper, LambdaTaper, CoreDiamOutput);
end

if fiber_type == 'Spiral'
    avgDistCores = LambdaTaper * 1.68; % empirical value from average distance of cores
else 
    avgDistCores = LambdaTaper;
end

%% Generate the fiber distal end - "EXPERIMENTAL VAIRABILITY"
% Polarisation randomisation: each core has a random linear polarisation
% angle theta_i uniform in [0, pi]. With a linear polariser in place,
% the transmitted amplitude is A_i = |cos(theta_i)|.
% Set SimulatePolarisation = 0 to disable (all cores at full amplitude).
SimulatePolarization = 0;
if SimulatePolarization
    theta_pol = pi * rand(n_cores, 1);        % random polarization angles
    polWeights = abs(cos(theta_pol));          % amplitude weights in [0,1]
    fprintf('Polarization randomisation enabled. Mean amplitude weight: %.3f (expected 2/pi = %.3f).\n', ...
        mean(polWeights), 2/pi);
else
    polWeights = ones(n_cores, 1);
    fprintf('Polarisation randomisation disabled.\n');
end
realAmpError = 0; % simulate not equal amplitudes in the core, between 0 and 1
if realAmpError ~= 0
    fprintf('Including amplitude errors of %.1f in the beamlets generation.\n', realAmpError);
end
GaussianAmpDecaySigma = 0; % simulate gaussian envelope for overfill of SLM
if GaussianAmpDecaySigma ~= 0
    fprintf('Including gaussian amplitude decay with sigma = %.1f um in the beamlets generation.\n', GaussianAmpDecaySigma);
end
fiber = FiberClassGenerator(fiber_type, n_cores, LambdaTaper);

%%%%%%%%%%%%%%%
%%% V & MFD %%%
%%%%%%%%%%%%%%%
EstimateParameters = 1;
if EstimateParameters == 1
    k0 = 2 * pi / lambda;
    if StepIndex == 1
        V = k0 * r_core * sqrt(n1^2 - n2^2); % STEP INDEX
        w0 = r_core * (0.65 + 1.619/(V^(3/2)) + 2.879 / (V^6) );
        MFD = w0 * 2;
        msg = sprintf('Estimaging V, MFD and NA using Marcuse formula for step index cores:\n');
        coreTypeString = 'StepIndexCore';
    else % GRIN fiber
        g = 2; % = 2 parabolic profile
        Delta = ( n1^2 - n2^2 ) / (2 * n1^2);
        V = k0 * n1 * r_core * sqrt(2*Delta); % GRADIENT INDEX
        msg = sprintf('V = %.2f', V);
        A = sqrt( 2/5 * ( 1 + 4 * (2/g)^(5/6) ) );
        B = exp(0.298/g) - 1 + 1.478 * (1 - exp(-0.077*g));
        C = 3.76 + exp(4.19 / g^(0.418) );
        w0 = r_core * ( A / (V^(2/(g+2)) ) + B / (V^(3/2)) + C / (V^6) );
        MFD  = 2 * w0;
        msg = sprintf('Estimaging V, MFD and NA using empirical formula for gradient index (GRIN) cores with g = %d:\n', g);
        coreTypeString = 'GradientIndexCore';
        % Ref: Gaussian approx. of fund. modes of graded-index fibers,
        % D. Marcuse 1977, JOSA
    end
else % manual MFD
    V = 3;
    MFD = 2.44 * 2;
    w0 = MFD / 2;
    msg = sprintf('MFD = %.2f microns manually selected.\n', MFD);
    coreTypeString = strcat('ManualParameters_MFD',num2str(MFD*10,'%d'));
end
NA_estimated = lambda / (pi * w0) ; % estimated from divergence of gaussian mode
z_r = (pi/lambda) *w0^2;
fill_factor = w0 / (avgDistCores/2);

disp(msg_setup);
disp(msg_taper);
fprintf('\nAverage distance between neighbour cores at output facet d_{avg} = %.2f um.\n', avgDistCores)
disp(msg);
fprintf('V = %.2f\nw0 = %.2f microns\nNA = %.2f\nz_r = %.2f um\n', V, w0, NA_estimated, z_r);
fprintf('\nThe focusing efficiency is affected by the fill factor at the output:\nw0/d_{avg} = %.2f\n', fill_factor);

%% Estimate minimum distance for all the cores to interfere
Dcores = zeros(n_cores, n_cores); % Initialize distance matrix
for i = 1:n_cores
    for j = 1:n_cores
        if i ~= j
            Dcores(i, j) = sqrt((fiber.posarray(i,1) - fiber.posarray(j,1))^2 + (fiber.posarray(i,2) - fiber.posarray(j,2))^2);
        else
            Dcores(i, j) = 0; % Set self-distance to infinity
        end
    end
end
D = max(Dcores,[],"all");
zmin = (D/2) * cot(asin(NA_estimated));
fprintf('Lambda_{out} = %.1f um, the max distance between cores is %.1f um.\n', LambdaTaper, D);
fprintf('Considering divergence NA = %.2f, the minimum distance for all the cores to interfere is is %d microns.\n', NA_estimated, uint16(zmin));


%% ===== FoV figures: model NA and (optionally) measured NA =====
pos    = fiber.posarray;          % um, output facet
amp    = polWeights(:);           % per-core amplitude weights
n_med  = 1;                       % 1 = air, 1.33 = water/tissue
z_show = z_focus;                  % plane shown in panel (c)
[~, idx] = max(Dcores(:));
[i1, i2] = ind2sub(size(Dcores), idx);

% --- Run 1: NA from the mode model (lambda / (pi*w0)) ---
resModel = makeFoVFigure(pos, D, i1, i2, lambda, NA_estimated, n_med, amp, z_show, ...
    sprintf('Model NA = %.3f (1/e^2), n_{med} = %.2f', NA_estimated, n_med));

% --- Run 2: experimentally measured NA ---
UseMeasuredNA = 1;
NA_meas_raw   = 0.147;             % value as measured, at the OUTPUT (tapered) facet
NA_meas_def   = '1/e2';           % '1/e2', '5%' or 'FWHM': definition used in the measurement
if UseMeasuredNA
    switch NA_meas_def
        case '1/e2', f = 1;
        case '5%',   f = sqrt(log(20)/2);   % 1.224
        case 'FWHM', f = sqrt(log(2)/2);    % 0.589
        otherwise, error('Unknown NA definition.');
    end
    NA_meas = NA_meas_raw / f;              % 1/e^2-equivalent NA (small-angle)
    resMeas = makeFoVFigure(pos, D, i1, i2, lambda, NA_meas, n_med, amp, z_show, ...
        sprintf('Measured NA = %.3f (%s) -> %.3f (1/e^2), n_{med} = %.2f', ...
        NA_meas_raw, NA_meas_def, NA_meas, n_med));

    fprintf('\n%-10s %8s %10s %16s %16s\n', 'Case', 'NA', 'z_min(um)', 'FoV@z_min(um)', ...
        sprintf('FoV@%dum(um)', z_show));
    fprintf('%-10s %8.3f %10.0f %16.1f %16.1f\n', 'Model', NA_estimated, ...
        resModel.zmin, resModel.FoV_zmin, resModel.FoV_zshow);
    fprintf('%-10s %8.3f %10.0f %16.1f %16.1f\n', 'Measured', NA_meas, ...
        resMeas.zmin, resMeas.FoV_zmin, resMeas.FoV_zshow);
end

%% ===================== LOCAL FUNCTIONS (keep at end of file) =====================
function res = makeFoVFigure(pos, D, i1, i2, lambda, NA, n_med, amp, z_show, figTitle)
    % Gaussian beamlet model fixed by the 1/e^2 NA
    lam_m = lambda / n_med;
    w0    = lambda / (pi * NA);                 % effective waist from NA
    zR    = pi * w0^2 / lam_m;
    wz    = @(z) w0 * sqrt(1 + (z/zR).^2);
    tanT  = tan(asin(NA / n_med));              % cone half-angle in the medium
    zmin  = (D/2) / tanT;
    n_cores = size(pos, 1);

    % Geometry for the side view
    c_pair = (pos(i1,:) + pos(i2,:)) / 2;
    u      = (pos(i2,:) - pos(i1,:)) / D;
    s      = (pos - c_pair) * u.';
    c_cent = mean(pos, 1);

    % Grid and envelope map (separable Gaussian sum)
    % zTop = 2*zmin;
    zTop = 600;
    if z_show > zTop
        zTop = 1.1*z_show;   % extend only if the selected plane lies beyond 2*z_min
        warning('z_show = %.0f um > 2*z_min; range extended to %.0f um.', z_show, zTop);
    end
    zz   = linspace(0.9*zmin, zTop, 80);
    L    = D/2 + 2.5*wz(zTop);
    dx   = 0.5;
    xv   = c_cent(1) + (-L:dx:L);
    yv   = c_cent(2) + (-L:dx:L);
    Imap = @(z) local_norm( ( exp(-(yv - pos(:,2)).^2 / wz(z)^2).' * ...
                              (amp .* exp(-(xv - pos(:,1)).^2 / wz(z)^2)) ).^2 );

    % Sweep FoV vs z
    hasIPT = ~isempty(which('bwconncomp'));
    FWHM_eq = zeros(size(zz)); FWHM_x = FWHM_eq; FWHM_y = FWHM_eq; nLobes = nan(size(zz));
    for m = 1:numel(zz)
        In = Imap(zz(m));
        [FWHM_eq(m), FWHM_x(m), FWHM_y(m)] = fwhmStats(In, dx);
        if hasIPT
            CC = bwconncomp(In >= 0.5); nLobes(m) = CC.NumObjects;
        end
    end
    zValid = zz(find(nLobes == 1, 1));

    % Figure
    figure('Color','w','Position',[60 60 1800 800]);
    t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    % t = tiledlayout(2, 3);
    % title(t, figTitle);

    % (a) Side view: cones only
    nexttile(1, [2 1]); hold on
    for k = 1:n_cores
        patch([s(k), s(k)-zTop*tanT, s(k)+zTop*tanT], [0 zTop zTop], ...
              [0.5 0.5 0.5], 'FaceAlpha', 0.02, 'EdgeColor', 'none');
    end
    % cols = [0.85 0.2 0.2; 0.2 0.3 0.85];
    cols = [0.2 0.3 0.85; 0.85 0.2 0.2];
    kk = [i1 i2];
    for m = 1:2
        k = kk(m);
        patch([s(k), s(k)-zTop*tanT, s(k)+zTop*tanT], [0 zTop zTop], ...
              cols(m,:), 'FaceAlpha', 0.2, 'EdgeColor', cols(m,:));
    end
    plot(s, zeros(size(s)), 'k.', 'MarkerSize', 8);
    yline(zmin,   'k--', sprintf('z_{min}'), 'LabelHorizontalAlignment', 'left', 'FontSize', 14);
    yline(z_show, 'r--',  sprintf('%.0f \\mum', z_show),     'LabelHorizontalAlignment', 'left', 'FontSize', 14);
    set(gca, 'YDir', 'reverse'); axis tight; box on
    xlabel('Position along farthest-core axis (\mum)'); ylabel('z (\mum)');
    % title('(a) Beamlet cones (1/e^2)');

    set(gca,'fontsize',16)

    % (b) Envelope at z_min
    nexttile(2);
    In_min = Imap(zmin);
    % plotEnvMap(xv, yv, In_min, pos, sprintf('(b) Focus envelope at z_{min} = %.0f \\mum', zmin));
    plotEnvMap(xv, yv, In_min, pos);
    xlim([-75 75])
    ylim([-75 75])

    set(gca,'fontsize',16)

    % (c) Envelope at z_show
    nexttile(5);
    In_show = Imap(z_show);
    % plotEnvMap(xv, yv, In_show, pos, sprintf('(c) Focus envelope at z = %.0f \\mum', z_show));
    plotEnvMap(xv, yv, In_show, pos);

    xlim([-75 75])
    ylim([-75 75])

    set(gca,'fontsize',16)

    % (d) FoV vs z
    nexttile(3, [2 1]); hold on; box on; grid on
    allF = [FWHM_eq FWHM_x FWHM_y];
    yl   = [max(0, min(allF) - 10), max(allF) + 10];
    if ~isempty(zValid) && zValid > zz(1)
        patch([zz(1) zValid zValid zz(1)], [yl(1) yl(1) yl(2) yl(2)], ...
              [0.9 0.9 0.9], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    plot(zz, FWHM_eq, 'k', 'LineWidth', 1.5);
    % plot(zz, FWHM_x, 'b--', zz, FWHM_y, 'r--');
    xline(zmin,   'k--', 'z_{min}', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off', 'FontSize', 14);
    xline(z_show, 'r--',  sprintf('%.0f \\mum', z_show), 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off', 'FontSize', 14);
    ylim(yl); xlim([zz(1) zz(end)]);
    xlabel('z (\mum)'); ylabel('FoV estimate (\mum)');
    % legend('Area-equivalent', 'x profile', 'y profile', 'Location', 'northwest');
    % title('(d) FoV vs distance');

    set(gca,'fontsize',16)
    % set(gcf(), 'DefaultAxesFontSize', 18)

    % Output
    res.zmin      = zmin;
    res.FoV_zmin  = fwhmStats(In_min,  dx);
    res.FoV_zshow = fwhmStats(In_show, dx);
    res.zz = zz; res.FWHM_eq = FWHM_eq; res.FWHM_x = FWHM_x; res.FWHM_y = FWHM_y;
end

function plotEnvMap(xv, yv, In, pos, varargin)
    if nargin > 4
        ttl = varargin{1};
    else
        ttl = 0;
    end
    imagesc(xv, yv, In); axis image; set(gca, 'YDir', 'normal'); colormap gray; hold on
    contour(xv, yv, In, [0.5 0.5], 'w', 'LineWidth', 1.5);
    plot(pos(:,1), pos(:,2), 'r.', 'MarkerSize', 7);
    clim([0 1]);                                  % caxis([0 1]) for releases before R2022a
    cb = colorbar; cb.Label.String = 'Relative focusing efficiency';
    xlabel('x (\mum)'); ylabel('y (\mum)');
    if ttl ~= 0
        title(ttl);
    end
end

function [d_eq, d_x, d_y] = fwhmStats(In, dx)
    d_eq = 2*sqrt(nnz(In >= 0.5)*dx^2/pi);
    [~, ip] = max(In(:)); [iy, ix] = ind2sub(size(In), ip);
    d_x = nnz(In(iy,:) >= 0.5) * dx;
    d_y = nnz(In(:,ix) >= 0.5) * dx;
end

function In = local_norm(I)
    In = I / max(I(:));
end