function Masks = phase_mask_mcf(basevecX, poscore, phi2, Nx, Ny, rot_angle)
% Generate the phase mask used for the MCF metasurface simulations.


Ncore = size(basevecX,1);
Nmeas = size(basevecX,2);

rot = [cosd(rot_angle) -sind(rot_angle); sind(rot_angle) cosd(rot_angle)];
poscore = round(poscore * rot);

[X,Y] = meshgrid(-Nx/2:Nx/2-1, -Ny/2:Ny/2-1);
[V,C] = voronoin(poscore);
IN = cell(1,Ncore);

for k = 1:Ncore
    IN{k} = find(inpolygon(X,Y,V(C{k},1),V(C{k},2)));
end

Masks = zeros(Ny,Nx,Nmeas);
for m = 1:Nmeas
    mask = zeros(Ny,Nx);
    for k = 1:Ncore
        idx = IN{k};
        mask(idx) = 2*pi*mod( ...
            angle(basevecX(k,m))/(2*pi) + ...
            phi2/(2*pi) * ((X(idx)-poscore(k,1)).^2 + (Y(idx)-poscore(k,2)).^2), 1);
    end
    Masks(:,:,m) = mask;
end
end
