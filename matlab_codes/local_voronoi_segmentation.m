function [IN, SegSizes] = local_voronoi_segmentation(poscore, Nx, Ny)
% Voronoi segmentation around the MCF core positions.

Ncore = size(poscore, 1);

    % Centered pixel grid
    [X, Y] = meshgrid(-Nx/2:Nx/2-1, -Ny/2:Ny/2-1);

    % Voronoi diagram
    [V, C] = voronoin(poscore);

    IN = cell(1, Ncore);
    SegSizes = zeros(Ncore, 1);

    for k = 1:Ncore

        vk = C{k};

        % Infinite Voronoi region or invalid region
        if any(vk == 1)
            IN{k} = [];
            SegSizes(k) = 0;
            continue;
        end

        tmp = V(vk, :);

        % Check invalid vertices
        if any(isinf(tmp(:))) || any(isnan(tmp(:)))
            IN{k} = [];
            SegSizes(k) = 0;
            continue;
        end

        % Pixels inside Voronoi polygon
        idx = find(inpolygon(X, Y, tmp(:,1), tmp(:,2)));

        IN{k} = idx;
        SegSizes(k) = numel(idx);

    end

end