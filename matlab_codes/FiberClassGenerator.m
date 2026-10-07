function fiber = FiberClassGenerator(fiber_type, n_cores, LambdaFiber, varargin)
    lambda = 0.920;
    fiber.type = fiber_type;
    A = ones(n_cores,1);
    relerrA = 0;         % VARARGIN(1), between 0 and 1
                         % relative error of core intensity
                         % can be used to model random polarization
                         % and error in coupling

    SigmaGaussDecay = 0; % sent as 7th input argument if you want 
                         %  exponential decay of the intensity as
                         % the cores go more far from center

    switch fiber_type
        case  {'reg','periodic', 'regular'}
        if n_cores > 1
            n_Rings = ( 1 + sqrt(1 + (4/3)*n_cores) ) / 2;
        else
            n_Rings = 1;
        end
        n_Rings = round(n_Rings);
        fprintf('Periodic fiber with %d rings\n', n_Rings);
        Mask = genMask('reg',n_Rings,LambdaFiber,0,0,0);
        n_cores = Mask.Nseg;
        posarray = genpos(Mask);
        fiber.posarray = posarray(1:n_cores,:);
        clear Mask
        otherwise
        fiber.posarray = fermatspiral(LambdaFiber, n_cores, 0, 0);
    end

    if nargin > 3
        relerrA = varargin{2};
        if relerrA > 1
            disp('Relative error can not be bigger than 1.')
            relerrA = 1;
        end
        if relerrA < 0
            disp('Relative error can not be smaller than 0.')
            relerrA = 0;
        end

        if nargin > 4
            SigmaGaussDecay = varargin{3};
            if SigmaGaussDecay > 0
                msg = sprintf('You selected an exponential decay of intensity with sigma %d microns.', SigmaGaussDecay);
                A = A .* exp(- (fiber.posarray(:,1).^2 + fiber.posarray(:,2).^2) / (2*SigmaGaussDecay^2) );
            else
                msg = sprintf('The sigma of intensityintensity decay needs to positive. Not applied.');
            end
            disp(msg);
        end
    end
    
    errA = relerrA * rand(size(A));
    A(2:end) = A(2:end) - errA(2:end); % first core is always the reference
    fiber.A = A;
end