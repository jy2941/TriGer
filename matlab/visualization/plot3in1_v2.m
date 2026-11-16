function plot3in1_v2(X, Y, XY, result, line, bar)
% Second version, put extra two below and on the left
% add p.val and cluster sepearate
% X-Y axis are reversed from imagesc function, this should be modified
% now assume m > n to simplify, I don't know what happens if n > m.......

    result = result(:, [2 1]);
    select_X = cat(2, result{:,1});
    select_Y = cat(2, result{:,2});
    % Step to remove duplicate part. 
    [~, idx] = unique(flip(select_X), 'stable');
    select_X = flip(select_X(end - idx + 1));
    [~, idx] = unique(flip(select_Y), 'stable');
    select_Y = flip(select_Y(end - idx + 1));

    C = Y(select_Y, select_Y);
    B = X(select_X, select_X);
    A = XY(select_X, select_Y);

    [m,n] = size(A);
    ratio = m/n/1.2; % 1.2 is kind of pre-set width-length ratio

    figure;
    
    %% --- Left: Diamonds for lower triangle
    ax1 = axes;

    % Count number of diamonds: sum of (n - i + 1) for i=1..n
    num = n*(n+1)/2;
    % Preallocate vertex arrays (each diamond has 4 vertices)
    plotX = zeros(4, num);
    plotY = zeros(4, num);
    plotC = zeros(1, num);
    
    vx = [0 1 0 -1];
    vy = [1 0 -1 0];
    
    count = 0;
    for i = 1:n
        for j = i:n
            count = count + 1;
            x0 = -i + 1;
            y0 = 2*n - j*2 + i;
    
            plotX(:,count) = x0 + vx;
            plotY(:,count) = y0 + vy;
            plotC(count) = C(j,j-i+1);
        end
    end
    
    % Build Faces matrix (each diamond = 4 vertices)
    Faces = reshape(1:4*num, 4, num)';
    
    % Combine all vertices into one big array (columns = vertices)
    Vertices = [plotX(:) plotY(:)];
    
    % Plot all diamonds in one patch call
    patch('Faces', Faces, 'Vertices', Vertices, ...
          'FaceVertexCData', plotC', 'FaceColor', 'flat', 'EdgeColor', 'none');
    axis equal tight off;

    %% --- Bottom: Mirrored diamonds for lower triangle ---
    
    num = m*(m+1)/2;
    
    % Preallocate vertex arrays (4 vertices per diamond)
    plotX = zeros(4, num);
    plotY = zeros(4, num);
    plotC = zeros(1, num);

    vx = [0 1/ratio 0 -1/ratio];
    vy = [1/ratio 0 -1/ratio 0];

    count = 0;
    for i = 1:m
        for j = i:m
            count = count + 1;
            x0 = (j*2 - i)/ratio;
            y0 = (-i + 1)/ratio;
    
            plotX(:,count) = x0 + vx;
            plotY(:,count) = y0 + vy;
            plotC(count) = B(j,j-i+1);
        end
    end
    
    % Build Faces matrix (each diamond = 4 vertices)
    Faces = reshape(1:4*num, 4, num)';
    
    % Combine all vertices into one big array (columns = vertices)
    Vertices = [plotX(:) plotY(:)];
    
    % Plot all diamonds in one patch call
    patch('Faces', Faces, 'Vertices', Vertices, ...
          'FaceVertexCData', plotC', 'FaceColor', 'flat', 'EdgeColor', 'none');

    %% --- Middle Square
    ax2 = axes;

    vxs = [0 2/ratio 2/ratio 0];
    vys = [0 0 2 2];
    
    num = m * n;  % total number of patches   
    % Preallocate
    plotX = zeros(4, num);
    plotY = zeros(4, num);
    plotC = zeros(1, num);
    
    count = 0;
    for i = 1:m
        for j = 1:n
            count = count + 1;
            x0 = (i - 1) * 2/ratio;
            y0 = (j - 1) * 2;
    
            plotX(:,count) = x0 + vxs;
            plotY(:,count) = y0 + vys;
            plotC(count) = A(i,n+1-j);
        end
    end
    
    % Faces and Vertices for one patch call
    Faces = reshape(1:4*num, 4, num)';
    Vertices = [plotX(:) plotY(:)];
    
    % Plot all patches at once
    patch('Faces', Faces, 'Vertices', Vertices, ...
          'FaceVertexCData', plotC', 'FaceColor', 'flat', 'EdgeColor', 'none');           % Set colormap
    axis equal tight off;
   %% Link axes
    linkaxes([ax1, ax2]); 
    caxis(ax1, [-1 1]); % colorbar for each plot
    if max(A(:)) > 1
        caxis(ax2, [0 25]);
    else
        caxis(ax2, [-1 1]);
    end    
    colormap(ax1, jet);
    colormap(ax2, jet);

    if nargin == 6
        cb = colorbar(ax2, 'eastoutside');
        cb.Position = [0.85 0.4 0.013 0.45];
    end

    %% Add line to separate blocks
    hold on 
    if line
        L=cellfun(@length, result); 
        L(:, 1) = L(:, 1) * 2 / ratio;
        L(:, 2) = L(:, 2) * 2;
        i=1;x=0;tmp=sum(L);y=tmp(2);
        while L(i, 1) ~= 0
           rectangle('Position', [x, y-L(i,2), L(i,1), L(i, 2)], 'EdgeColor', 'k', 'LineWidth', 2);
           y=y-L(i,2);
           x=x+L(i,1);
           i = i+1;
        end
    end
    hold off
end    
