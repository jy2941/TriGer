clear;
addpath('../matlab/methods')
rng(24);

% Define the size of the matrix, and generate element matrix
% X: 450, Y:550
size = 1000;
R = eye(size);
m = 450;
n = 550;

% 40 correlated X, 60 correlated Y. 

for i=1:40
    for j=451:510
        R(i,j) = 0.3; R(j,i) = 0.3; 
    end 
end 

for i=1:40
    for j=1:40
        if j>i
            R(i,j) = 0.5; R(j,i) = 0.5;
        end
    end    
end

for i=451:510
    for j=451:510
        if j>i
            R(i,j) = 0.5; R(j,i) = 0.5;
        end
    end    
end

chol(R);

% Second Part, 40 correlated X and 30 extra correlated Y

for i=1:40;
    for j=511:540
        R(i,j) = 0.2; R(j,i) = 0.2;
    end   
end 


for i=511:540
    for j=511:540
        if j>i
            R(i,j) = 0.5; R(j,i) = 0.5;
        end
    end    
end

R = (R+R')/2;
figure;imagesc(R);colorbar;colormap(jet(256)); caxis([-0.2, 1]); ax=gca;ax.FontSize=18;
chol(R);

clearvars -except R m n size
save('matrix_1.mat');




%% One time simulation with TGI
mu = zeros(1, size); % Mean vector of zeros
data = mvnrnd(mu, R, 2000);
[sample_R, sample_P] = corr(data, 'Rows', 'pairwise');
figure;imagesc(sample_R);colorbar;colormap(jet(256)); caxis([-0.2, 1]);ax=gca;ax.FontSize=18;

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
figure;imagesc(shuffle_R);colorbar;colormap(jet(256)); caxis([-0.2, 1]);ax=gca;ax.FontSize=18;
hold on;
line([0 1000], [450 450], 'Color', 'r', 'LineWidth', 2);
line([450 450], [0 1000], 'Color', 'r', 'LineWidth', 2);
hold off;

% test it can recover from org_idx
figure;imagesc(shuffle_P(org_idx, org_idx));colorbar;colormap(jet(256)); caxis([0, 4]);ax=gca;ax.FontSize=18;
shuffle_R = abs(shuffle_R);

% Start TGI method
figure;imagesc(shuffle_R);colorbar;colormap(jet(256)); caxis([-0.2, 1]); ax=gca;ax.FontSize=18;
result = greedy_3step(shuffle_R, shuffle_R, m, n, 1, 1.8, 1.8, [10, 50], [10, 50])

