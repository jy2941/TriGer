function result = sort_result(Wp, result)
% apply for X-Y matrix only, wait for update
% now apply for p*p matrix!
    m = size(Wp, 1);
    n = size(Wp, 2);   
    i = 1;
    if m == n
        while ~isempty(result{i, 1})
            sub_Wp = Wp(result{i, 1}, result{i, 1});
            R = sum(sub_Wp,2);
            [~, order] = sort(R, 'descend');
            result{i, 1} = result{i, 1}(order);
            i = i+1;
        end    
    else
        while ~isempty(result{i, 1})
            sub_Wp = Wp(result{i, 1}, result{i, 2});
            C = sum(sub_Wp,1); 
            R = sum(sub_Wp,2);
            [~, order] = sort(R, 'descend');
            result{i, 1} = result{i, 1}(order);
            [~, order] = sort(C, 'descend');
            result{i, 2} = result{i, 2}(order);
            i = i+1;
        end
        % For reorder the cluster, comment it now. 
        % nRows = size(result, 1);
        % vec = zeros(1, nRows);
        % for i = 1:nRows
        %    if numel(result{i,1}) > 0
        %        sub = Wp(result{i,1}, result{i,2});
        %        vec(i) = mean(sub(:));
        %    end    
        % end
        % Get the sorted order of mean vector
        % [~, order] = sort(vec, 'desc');
        % Reorder result using the order
        % result = result(order, :);
    end
end


