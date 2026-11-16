function result = sort_result_XandY(X, Y, result)
% apply for X-Y matrix only, wait for update
% sort the result, according to X&Y covariance matrix, not X-Y matrix! 
    i = 1;
    while ~isempty(result{i, 1})
        sub_X = X(result{i, 1}, result{i, 1});
        C = sum(sub_X,1); 
        [~, order] = sort(C, 'descend');
        result{i, 1} = result{i, 1}(order);
        sub_Y = Y(result{i, 2}, result{i, 2});
        R = sum(sub_Y,1);
        [~, order] = sort(R, 'descend');
        result{i, 2} = result{i, 2}(order);
        i = i+1;
    end
end