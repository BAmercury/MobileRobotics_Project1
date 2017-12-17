function plotEnv(map, blocks)
    % Draw the occupancy grid
    map.show();
    % Overlay the blocks
    hold on;
    block_plot(gcf, blocks(:, 1:2)', blocks(:, 3:4)', [255 0 0]');
    hold off;
end

function block_plot(h, p1, p2, c)
% Plot a list blocks on the given figure handle h,
% p1, xmin ymin zmin
% p2, xmax ymax zmax
% c, color in r g b (0-255)

    v = [0 0; ...
         1 0; ...
         1 1; ...
         0 1];
    f = [1 2 3 4];
    c = c / 255;
    
    if size(c,2) == 1
        c = repmat(c, 1, size(p1,2));
    end

    d   = p2 - p1;
    N   = size(p1,2);
    vl  = zeros(4*N,2);
    fl  = zeros(N,4);
    fcl = zeros(N,3);
    for k = 1:size(p1,2)
        vs = v;
        vs(:,1) = vs(:,1) * d(1,k);
        vs(:,2) = vs(:,2) * d(2,k);
        vs(:,1) = vs(:,1) + p1(1,k);
        vs(:,2) = vs(:,2) + p1(2,k);
        vl((k-1)*4+1:k*4,:) = vs;
        fl(k,:) = f + (k-1)*4;
        fcl(k,:) = c(:,k)';
    end

    figure(h);
    hold on;
    patch('Vertices', vl, 'Faces', fl, 'FaceColor', 'flat', ...
        'FaceVertexCData', fcl, 'EdgeAlpha', 0, 'FaceAlpha', 0.6);
    hold off;
end