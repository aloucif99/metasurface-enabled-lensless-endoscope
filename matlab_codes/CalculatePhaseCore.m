function phi = CalculatePhaseCore(Xc,Yc, target_pos, lambda, varargin)
    rel_phase_err = 0; 

    % geometrical factors
    phi1x = (2*pi/lambda) * target_pos(1) / target_pos(3);
    phi1y = (2*pi/lambda) * target_pos(2) / target_pos(3);
    phi2 = - (2*pi/lambda) / target_pos(3);

    phi = phi1x * Xc + phi1y * Yc + (1/2) * phi2 * (Xc^2 + Yc^2);
    
    if nargin > 4
        rel_phase_err = varargin{1}; % random field phase error relative to pi (1 = pi error, 0 = no error)
        phase_err = rel_phase_err * pi * (1-2*rand(1));
        phi = phi + phase_err;
    end

end
