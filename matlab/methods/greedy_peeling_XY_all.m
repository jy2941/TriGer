function result = greedy_peeling_XY_all(XY, lambda, max_round)
% Simpler greedy algorithm based only on X-Y matrix
% input: XY matrix, and lambda
% extract ALL submatrix
% Default round is 5

    if  nargin < 3
        max_round = 5;
    end

    result = cell(10, 2);
    i = 1;

    % Original Density 
    m = size(XY, 1);
    n = size(XY, 2);
    org_density = sum(sum(XY)) / (m*n)^(lambda/2); 

    remain_X = 1:m;
    remain_Y = 1:n;
    Wp = XY;

    while true
    
        disp(['Round ', num2str(i)])
        disp('******************')

        [W_density, select_X, select_Y]=greedy_peeling_XY_one(Wp, lambda);

        if W_density > org_density & i < max_round
            result{i, 1} = remain_X(select_X);
            result{i, 2} = remain_Y(select_Y);
            all_select_X = cat(2, result{:,1});
            all_select_Y = cat(2, result{:,2});
            remain_X = setdiff(1:m, all_select_X);
            remain_Y = setdiff(1:n, all_select_Y);
            Wp = XY(remain_X, remain_Y);
            i = i+1;
        else 
            break
        end    
    end 
end