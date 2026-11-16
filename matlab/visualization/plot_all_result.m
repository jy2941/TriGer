function plot_all_result(Wp, result, full, line, crange)
% plot all regions from greedy peeling
% Input matrix should be X (symmetric) or X-Y
% full refers whole matrix or only selected region
% line refers whether to add line (region out of whole matrix / divide region)
% also input color range

    m = size(Wp, 1);
    n = size(Wp, 2);

    if m == n % Manually assign second row in result, make it consistent
        result(:, 2) = result(:, 1);
    end    

    select_X = cat(2, result{:,1});
    select_Y = cat(2, result{:,2});

    if full % extract matirx according to whether full or region
        remain_X = setdiff(1:m, select_X);
        remain_Y = setdiff(1:n, select_Y);
        order_X = [select_X, remain_X];
        order_Y = [select_Y, remain_Y];
        Wp_reordered = Wp(order_X, order_Y);
    else 
        Wp_reordered = Wp(select_X, select_Y);
    end     


    if nargin < 5 % wheter color axis is provided
        crange = [min(Wp_reordered(:)) max(Wp_reordered(:))];
    end
       
    figure;imagesc(Wp_reordered);colorbar;colormap(jet);caxis(crange);
    if line
        hold on
        if full
            rectangle('Position', [0, 0, ...
                length(select_Y) + 0.5, length(select_X) + 0.5], 'EdgeColor', 'r', 'LineWidth', 2);
        else 
            L = cellfun(@length, result);
            i=1;x=0;y=0;
            L(i,1)=L(i,1)+0.5;L(i,2)=L(i,2)+0.5;
            while L(i, 1) ~= 0
               rectangle('Position', [x, y, L(i,2), L(i,1)], 'EdgeColor', 'r', 'LineWidth', 2);
               y=y+L(i,1);
               x=x+L(i,2);
               i = i+1;
            end
        end
        hold off
    end    
end 