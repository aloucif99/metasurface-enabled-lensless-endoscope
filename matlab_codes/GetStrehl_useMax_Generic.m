function [Strehl, varargout] = GetStrehl_useMax_Generic(Img, varargin)
    posMax = zeros(2,1);
    [~, ind] = max(Img,[],'all');
    [posMax(1), posMax(2)] = ind2sub(size(Img),ind);
%     [posMax(2), posMax(1)] = ind2sub(size(Img),ind);
    range_fit = 40;
    CheckWhatsGoinOn = 0;
    if nargin > 1
        CheckWhatsGoinOn = varargin{1};
        if nargin > 2
                range_fit = varargin{2};
        end
    end

    range_array = -range_fit/2:range_fit/2;
    Xprofile = Img(posMax(1)+range_array,posMax(2));
    Yprofile = Img(posMax(1),posMax(2)+range_array);

    fx = fit(range_array',Xprofile-min(Xprofile),'gauss1');
    fy = fit(range_array',Yprofile'-min(Yprofile),'gauss1');
    
    sigma_x = fx.c1 / sqrt(2);
    sigma_y = fy.c1 / sqrt(2);
    sigma = mean([sigma_x sigma_y]);
    FWHM = 2 * sqrt(2*log(2)) * sigma;

    %% create binary mask to identify if a pixel is inside the focus
    [Y,X] = meshgrid(1:size(Img,1),1:size(Img,2));
    Nsigma = 3;
    RadiusFocus = Nsigma*sigma;
    DistanceFromFocus = sqrt( ( X-posMax(1) ).^2 + (Y - posMax(2)).^2 );
    FocusMap = zeros(size(Img));
%     FocusMap(DistanceFromFocus <= RadiusFocus) = 1;

    for i = 1:size(Img,1)
        for j = 1:size(Img,2)
            if (posMax(1) - X(i,j) )^2 +  (posMax(2) - Y(i,j) )^2 <= (RadiusFocus)^2
                FocusMap(i,j) = 1;
            end
        end
    end
 
    TotalPower = sum(Img,"all");
    FocusPower = sum(Img.*FocusMap,"all");
    
    Strehl = FocusPower/TotalPower;

    if CheckWhatsGoinOn
        figure(5)
        subplot(2,1,1)
        plot(posMax(1)+range_array,Img(posMax(1)+range_array,posMax(2)))
        xlabel('X [pixels]')
        ylabel('Intensity')
        subplot(2,1,2)
        plot(posMax(2)+range_array,Img(posMax(1),posMax(2)+range_array))
        xlabel('Y [pixels]')
        ylabel('Intensity')
        
        figure(6)
        surf(Y, X, Img .* FocusMap)
        view(2), shading interp, axis square, colormap("jet")
        title(strcat('Image * FocusMask with  ',num2str(Nsigma,'%d'), ' \sigma from gaussian fit'))
        pbaspect([1 1 1])
    end

    if nargout > 1
    varargout{1} = posMax;
        if nargout > 2
            varargout{2} = FWHM;
            if nargout > 3
                varargout{3} = FocusPower;
            end
        end
    end
end