function plot3in1_v1(X, Y, XY, result, line, step, bar)
% first version, extra two are top and on the right
% now assume m > n to simplify, don't ask me what happens if n > m.......
% result is a p*2 cell, each row records X&Y index
% line is true or false, whether add cluster line
% step could be ignored, step size for X&Y tick, if not assigned, 10&100 by default
% bar is colorbar, could be ignored
    
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

    %% --- Right: Diamonds for lower triangle
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
            x0 = i - 1 + m*2/ratio;
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
    %% --- Top: Mirrored diamonds for lower triangle ---
    
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
            y0 = (i - 1)/ratio + n*2;
    
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
    axis equal tight off;

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
    axis equal tight;

    box on;
    ax2.Color = 'none';
    ax2.XRuler.Axle.LineStyle = 'none';  % hides x-axis line
    ax2.YRuler.Axle.LineStyle = 'none';  % hides y-axis line

    % --- Clean Tick Labels ---
    clean_ticks = @(maxVal, base) base:base:floor(maxVal/base)*base;
    
    if nargin < 6
        xs = 100; ys = 10;
    else
        xs = step(1);ys=step(2);
    end    
    xTickIdx = clean_ticks(m, xs);  % X: 20, 40, ..., 200
    yTickIdx = flip(clean_ticks(n, ys));  % Y: 10, 20, ..., 50

    xTicks = (xTickIdx - 1) * 2/ratio;
    yTicks = (n - yTickIdx) * 2 + 1;

    % Set ticks with reversed Y tick labels (descending)
    set(gca, 'XTick', xTicks, 'XTickLabel', string(xTickIdx));
    set(gca, 'YTick', yTicks, 'YTickLabel', string(yTickIdx));
    box off;

   %% Link axes
    linkaxes([ax1, ax2]); 
    caxis(ax1, [-1 1]); % colorbar for each plot
    if max(A(:)) > 1
        caxis(ax2, [0 10]); % Maximum value for middle plot, I would suggest modify by hand
    else
        caxis(ax2, [-1 1]);
    end    
    colormap(ax1, jet);
    colormap(ax2, jet);

    if nargin == 7
        cb = colorbar(ax2, 'eastoutside');
        cb.Position = [0.15 0.13 0.013 0.45];
    end

    %% Add line to separate blocks
    hold on 
    if line
        L=cellfun(@length, result); 
        L(:, 1) = L(:, 1) * 2 / ratio;
        L(:, 2) = L(:, 2) * 2;
        i=1;x=0;tmp=sum(L);y=tmp(2);
        while L(i, 1) ~= 0
           rectangle('Position', [x, y-L(i,2), L(i,1), L(i, 2)], 'EdgeColor', 'r', 'LineWidth', 2);
           y=y-L(i,2);
           x=x+L(i,1);
           i = i+1;
        end
    end
    hold off
end    
