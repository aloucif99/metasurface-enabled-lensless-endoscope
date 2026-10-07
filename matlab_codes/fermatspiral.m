function [posarray] = fermatspiral(d,n,theta,SetFirstToCenter)
%% Fermat's spiral
% returns the cartesian co-ordinates of the points describing a fermat's
% spiral centered at [0,0]
%
% d == distance between the spots
%
% n == number of spots
%
% Notation from Aperiodic Antenna Array for Secondary Lobe Suppression,
% IEEE PHOTONICS TECHNOLOGY LETTERS, VOL. 28, NO. 2, JANUARY 15, 2016
%
% $Sid Last edit :: 10/11/2016
% Luca :: added the theta angle and the set of the fermat to center

if nargin<1
    d = 20;
    n = 169;
    theta = 0;
    SetFirstToCenter = 0;
    Clockwise = 0;
else
    if nargin<3
        theta = 0;
        SetFirstToCenter = 0;
        Clockwise = 0;
    else
        if nargin<4
            SetFirstToCenter = 0;
            Clockwise = 0;
        end
    end
end

% d14 = sqrt(5-4*cos(3*pi*(3-sqrt(5)))); % normalization factor (I don't understand clearly why)
d14 = 1;
n_array = 1:n;

rho = d/d14*sqrt(n_array);
phi = ( n_array*pi*(3-sqrt(5)) );

[X,Y] = pol2cart(phi + theta,rho);


if SetFirstToCenter == 1
    posarray = [X-X(1); Y-Y(1)]';
else
    posarray = [X; Y]';
end

% scatter(X,Y)
% figure
end