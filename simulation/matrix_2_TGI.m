addpath('../matlab/methods')
load('matrix_1.mat')

nRuns  = 10;
TPR_all = zeros(nRuns, 1);
TNR_all = zeros(nRuns, 1);
sample_n = 500;

for run = 1:nRuns
    rng(run);
    %% Now we have the population level correlation matrix, generate sample correlation matrix 
    mu = zeros(1, size); % Mean vector of zeros
    data = mvnrnd(mu, R, sample_n);
    [sample_R, sample_P] = corr(data, 'Rows', 'pairwise');  
    %% shuffle the sample correlation matrix    
    shuffle_X_idx = randperm(450); % Generate random permutation of X
    [~, org_X_idx] = sort(shuffle_X_idx);
    shuffle_Y_idx = randperm(550) + 450; % Generate random permutation of Y
    [~, org_Y_idx] = sort(shuffle_Y_idx);
    org_idx = [org_X_idx, org_Y_idx+450];
    shuffle_idx = [shuffle_X_idx, shuffle_Y_idx];
    
    sample_R(1:size+1:end) = 0; 
    shuffle_R = sample_R(shuffle_idx, shuffle_idx); % Apply permutation to rows and columns
    shuffle_P = sample_P(shuffle_idx, shuffle_idx);
    shuffle_P = -log10(shuffle_P);
    shuffle_R = abs(shuffle_R);
    
    % Start TGI method
    result = greedy_3step(shuffle_R, shuffle_R, m, n, 1.2, 1.8, 1.8, [10, 50], [10, 50])
    
    % Map result back to original space
    detected_X = result{1,1};
    detected_Y = result{1,2} - m;  % remove m offset before indexing

    % True/false X (original space)
    true_X  = org_X_idx(1:40);
    false_X = org_X_idx(41:450);
    
    true_Y  = org_Y_idx(1:90);    
    false_Y = org_Y_idx(91:550);
    
    % TPR/TNR for X
    TP_X = numel(intersect(detected_X, true_X));
    TN_X = numel(setdiff(false_X, detected_X));
    TPR_X = TP_X / 40;
    TNR_X = TN_X / 410;
    
    % TPR/TNR for Y
    TP_Y = numel(intersect(detected_Y, true_Y));
    TN_Y = numel(setdiff(false_Y, detected_Y));
    TPR_Y = TP_Y / 90;
    TNR_Y = TN_Y / 460;
    
    % Average TPR and TNR across X and Y
    TPR_all(run) = (TPR_X + TPR_Y) / 2;
    TNR_all(run) = (TNR_X + TNR_Y) / 2;
end