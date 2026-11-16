function result = greedy_peeling_X_all(Wp, lambda)
    % Simpler greedy algorithm based only on X-X matrix
    % input: symmetric matrix, and lambda
    % extract ALL submatrix

    if ~issymmetric(Wp)
        error('Matrix Input must be symmetric.');
    end  

    result = cell(10, 2);
    i = 1;
    
    % Original Density 
    m = size(Wp, 1);
    org_density = sum(sum(Wp)) / (m^lambda); 
    X = Wp; % Original Matrix backup

    remain_X = 1:m; % Avaiable nodes

    while true
    
        disp(['Round ', num2str(i)])
        disp('******************')

        [W_density, select_X]=greedy_peeling_X_one(Wp, lambda);

        if W_density > org_density
            result{i, 1} = remain_X(select_X); % selected nodes, idx for original Wp
            result{i, 2} = W_density;
            all_select_X = cat(2, result{:,1});
            remain_X = setdiff(1:m, all_select_X);
            Wp = X(remain_X, remain_X);
            i = i+1;
        else 
            break
        end    
    end 
end