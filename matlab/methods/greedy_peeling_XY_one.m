function [W_density, select_X, select_Y] = greedy_peeling_XY_one(Wp, lambda)
    % Simpler greedy algorithm based only on X-Y matrix
    % extract only ONE submatrix
    % output is only the greedy matrix
    
    Wp_temp = single(Wp);
    W_density = 0;
    m = size(Wp_temp, 1);
    n = size(Wp_temp, 2);
    remove_X = [];
    remove_Y = [];
    remove_X_final = [];
    remove_Y_final = [];
    len_X = m;
    len_Y = n;
    C = sum(Wp_temp,1);  % column sum
    R = sum(Wp_temp,2);  % row sum
    for i = 1:(m+n)
        C(remove_Y) = inf;  % Exclude already removed indices for Y
        R(remove_X) = inf;  % Exclude already removed indices for X
        [dT, IndT] = min(C);   
        [dS, IndS] = min(R);  
        if dT == inf | dS == inf   % stop the loop if all X or all Y are excluded. 
            break;
        end
        
        % Update remove_idx
        c = len_X/len_Y;
        if c*dS <= dT
            remove_X(end+1) = IndS;
            C = C - Wp_temp(IndS, :);
            C(remove_Y) = 0;
            Wp_sum = sum(C);
            len_X = len_X - 1;
        else
            remove_Y(end+1) = IndT;
            R = R - Wp_temp(:, IndT);
            R(remove_X) = 0;
            Wp_sum = sum(R);
            len_Y = len_Y - 1;
        end
        
        score_temp = Wp_sum / (len_X*len_Y)^(lambda/2);
        if score_temp > W_density & score_temp ~= Inf
            W_density = score_temp;
            remove_X_final = remove_X;
            remove_Y_final = remove_Y;
        end    
        % disp(i)
    end
    
    % Find the index of the maximum score
    select_X = setdiff(1:m, remove_X_final);
    select_Y = setdiff(1:n, remove_Y_final);
end