function result = greedy_3step(P, R, m, n, lambda1, lambda2, lambda3, X_size, Y_size, min_community_size)
% Full 3-step greedy peeling.
% P: (m+n) x (m+n) association matrix (e.g. -log10 p); P(1:m, m+1:m+n) drives step 1
% R: (m+n) x (m+n) signed correlation matrix for within-omics refinement
% X_size, Y_size: [min max]. Blocks with min or fewer features are dropped;
% blocks with more than max features are refined by within-omics peeling.
% Refined communities smaller than min_community_size are dropped; if none is
% kept, that side of the block stays unrefined.
% result: p x 2 cell, {X indices, Y indices} (Y indices in 1..n)

    if nargin < 5, lambda1 = 1.5; end
    if nargin < 6, lambda2 = 1.5; end
    if nargin < 7, lambda3 = 1.5; end
    if nargin < 8, X_size = [10 200]; end
    if nargin < 9, Y_size = [10 200]; end
    if nargin < 10, min_community_size = 2; end

    blocks = greedy_peeling_XY_all(P(1:m, m+1:m+n), lambda1);
    result = cell(0, 2);
    i = 1;
    while i <= size(blocks, 1) && ~isempty(blocks{i, 1})
        xi = blocks{i, 1}; yi = blocks{i, 2};
        i = i + 1;
        if numel(xi) <= X_size(1) || numel(yi) <= Y_size(1)
            continue
        end
        if numel(xi) > X_size(2)
            xi = refine(xi, R(xi, xi), lambda2, min_community_size);
        end
        if numel(yi) > Y_size(2)
            yi = refine(yi, R(m + yi, m + yi), lambda3, min_community_size);
        end
        result(end + 1, :) = {xi, yi};
    end
end

function idx = refine(idx, sub, lambda, min_size)
    sub = (sub + sub') / 2;
    comm = greedy_peeling_X_all(sub, lambda);
    comm = comm(~cellfun(@isempty, comm(:, 1)), :);
    comm = comm(cellfun(@numel, comm(:, 1)) >= min_size, :);
    if isempty(comm)
        return
    end
    for k = 1:size(comm, 1)
        c = comm{k, 1};
        [~, order] = sort(sum(sub(c, c), 2), 'descend');
        comm{k, 1} = c(order);
    end
    idx = idx(cat(2, comm{:, 1}));
end
