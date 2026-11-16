function [W_density, select_X] = greedy_peeling_X_one(Wp, lambda)
    % Simpler greedy algorithm based only on X-X matrix
    % Input must be a symmetric matrix
    % extract only ONE submatrix
    % output is only the greedy matrix
    Wp_temp = single(Wp);

    if ~issymmetric(Wp_temp)
        error('Matrix A must be symmetric.');
    end    

    m = size(Wp_temp, 1);
    W_density = 0;
    remove_X = [];
    remove_X_final = [];
    C = sum(Wp_temp, 1);  % column sum
    for i = 1:m
        C(remove_X) = inf;  % Exclude already removed indices for X
        [~, IndT] = min(C);         
        % Update remove_idx
        remove_X(end+1) = IndT;
        C = C - Wp_temp(IndT, :);
        C(remove_X) = 0;
        Wp_sum = sum(C);        
        score_temp = Wp_sum / (m-i)^(lambda);
        if score_temp > W_density & score_temp ~= Inf
        % Density > prior density and not Inf
            W_density = score_temp;
            remove_X_final = remove_X;
        end    
        % disp(i)
    end
    select_X = setdiff(1:m, remove_X_final);
end