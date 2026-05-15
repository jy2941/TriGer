%% Dense Community Detection in Multi-omics Data (Dataset 2)
%
% This script reproduces all figures and results for Dataset 2.
% Input: CleanCorNonzero.mat — pre-computed correlation matrices
% (Pearson r and -log10 p-value for Treatment and Control groups).
% No raw expression data is used.
%
% *Analysis pipeline:*
%
% # Load correlation matrices (R_Treat, R_Control, P_Treat, P_Control)
% # Combine into a single R and P matrix (element-wise max by magnitude)
% # Build bipartite weight matrix Wp = P(X-genes, Y-genes), threshold p < 0.01
% # Run greedy peeling to detect dense bipartite communities (clusters)
% # Refine Cluster 1 by re-running peeling on each side independently
% # Split and reorder clusters, generate heatmaps and summary statistics

%% Setup
addpath('methods')
addpath('visualization')

set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultTextFontSize', 18);
set(0, 'DefaultAxesColor', 'white');
set(0, 'DefaultFigureColor', 'white');

%% Load Correlation Matrices
% CleanCorNonzero.mat contains four pre-computed matrices (single precision):
%   R_Treat, R_Control  — Pearson correlation coefficients
%   P_Treat, P_Control  — significance as -log10(p-value)
%   name                — cell array of gene/feature names (1 x 1428)
%
% Dimensions: 341 X-features (rows 1:341) and 1087 Y-features (rows 342:1428).
load('Real_data/CRC/CleanCorNonzero.mat');

m = 341;
n = 1087;

fprintf('Loaded CleanCorNonzero.mat\n');
fprintf('  X features : %d\n', m);
fprintf('  Y features : %d\n', n);
fprintf('  Matrix size: %d x %d\n', m+n, m+n);

%% Build Combined Correlation and Significance Matrices
% P = element-wise max of P_Treat and P_Control.
%   Keeps the strongest significance signal from either condition.
%
% R = element-wise max by absolute value between R_Treat and R_Control.
%   Where Control has a larger absolute correlation, it replaces Treat.
%   This ensures the combined R reflects the dominant correlation structure.
P = max(P_Treat, P_Control);

R = R_Treat;
mask = abs(R_Control) >= abs(R_Treat);
R(mask) = R_Control(mask);

%% Construct Bipartite Weight Matrix
% Wp = P(X-genes, Y-genes), entries below -log10(0.01) = 2 set to zero.
% This is the bipartite graph on which community detection runs.
Wp = P(1:m, (m+1):(m+n));
Wp(Wp < 2) = 0;

figure;
imagesc(Wp); colorbar; colormap(jet); caxis([0 15]);
title('Bipartite weight matrix W_p  (X \times Y)');
xlabel('Y gene index'); ylabel('X gene index');

%% Initial Community Detection: Greedy Peeling
% greedy_peeling_XY_all finds dense bipartite subgraphs by iteratively
% peeling away low-weight nodes. Threshold 1.42 controls minimum density.
%
% result is an Nx2 cell array:
%   result{i,1} = X member indices for cluster i
%   result{i,2} = Y member indices for cluster i  (relative to 1, not m+1)
result = greedy_peeling_XY_all(Wp, 1.42);

% Sort clusters by internal X and Y correlation structure.
result = sort_result_XandY(R_Treat(1:m, 1:m), R_Treat((m+1):(m+n), (m+1):(m+n)), result);
result([3 4], :) = result([4 3], :);
result = sort_result(P(1:m, (m+1):(m+n)), result);

n_clusters_init = sum(~cellfun(@isempty, result(:,1)));
fprintf('Initial clusters detected: %d\n', n_clusters_init);

%% Initial Result: Bipartite Weight Matrix Reordered by Clusters
plot_all_result(P(1:m, (m+1):(m+n)), result, true,  true, [0 15]);
plot_all_result(P(1:m, (m+1):(m+n)), result, false, true, [0 15]);

%% Refine Cluster 1: Sub-community Detection
% Cluster 1 is the largest community. Re-running peeling on its X and Y
% sides separately reveals finer sub-structure within each modality.
idx_X = cat(1, result{1,1});
idx_Y = m + cat(2, result{1,2});

X = R(idx_X, idx_X);
Y = R(idx_Y, idx_Y);

% Re-peel X side of Cluster 1
test = greedy_peeling_X_all(X, 1.3);
test = sort_result(X, test);
result{1,1} = idx_X(cat(2, test{:,1}));

% Re-peel Y side of Cluster 1
test = greedy_peeling_X_all(Y, 1.3);
test = sort_result(Y, test);
result{1,2} = idx_Y(cat(2, test{:,1})) - m;

%% Split Cluster 1 into Two Sub-clusters
% After refinement, Cluster 1 contains distinguishable sub-groups.
% Rows 18-25 of X (and rows 107-197 of Y) form a separate sub-cluster.
result = [result(1,:); cell(1,2); result(2:end,:)];
result{2,1} = result{1,1}(18:25);
result{1,1} = result{1,1}(1:17);
result{2,2} = result{1,2}(107:197);
result{1,2} = result{1,2}(1:106);

% Swap Y memberships between sub-clusters 1 and 2 for visual alignment
tmp         = result{1,2};
result{1,2} = result{2,2};
result{2,2} = tmp;

%% Reorder Remaining Clusters
result([3 4 5], :) = result([4 5 3], :);
result = sort_result_XandY(R_Treat(1:m, 1:m), R_Treat((m+1):(m+n), (m+1):(m+n)), result);

n_clusters = sum(~cellfun(@isempty, result(:,1)));
fprintf('Final clusters after refinement: %d\n', n_clusters);

%% 3-Panel Overview: Y-Y, X-X, and Y-X Correlation
% Shows all three correlation blocks reordered by the final cluster assignment.
%   Left  : Y intra-correlation
%   Middle: X intra-correlation
%   Right : Y-X cross-correlation
X_mat  = R_Treat((m+1):(m+n), (m+1):(m+n));
Y_mat  = R_Treat(1:m, 1:m);
XY_mat = R_Treat((m+1):(m+n), 1:m);
plot3in1_v1(X_mat, Y_mat, XY_mat, result, true, [50, 10], true);

%% Per-modality Heatmaps Reordered by Clusters
plot_all_result(R(1:m, (m+1):(m+n)), result,         false, true, [-1,  1]);
plot_all_result(R(1:m, 1:m),         result(:,1),     false, true, [-1,  1]);
plot_all_result(R((m+1):(m+n), (m+1):(m+n)), result(:,2), false, true, [-1, 1]);

%% Full Matrix Comparison: Combined, Control, Treat
% Visualise the same cluster reordering across three versions of R
% to compare how the structure looks in each condition.
for Rt = {R, R_Control, R_Treat}
    Rt = Rt{1};
    X_mat  = Rt((m+1):(m+n), (m+1):(m+n));
    Y_mat  = Rt(1:m, 1:m);
    XY_mat = Rt((m+1):(m+n), 1:m);
    plot3in1_v1(X_mat, Y_mat, XY_mat, result, true, [50, 10], true);
end
plot_all_result(R(1:m, (m+1):(m+n)), result, false, true, [-1, 1]);

%% Cluster Size and Treat vs Control Mean Correlation
% L(i,1) = mean r across all X-Y gene pairs in cluster i, Control group
% L(i,2) = mean r across all X-Y gene pairs in cluster i, Treat group
XY_Treat   = R_Treat(1:m, (m+1):(m+n));
XY_Control = R_Control(1:m, (m+1):(m+n));

L = [];
i = 1;
while ~isempty(result{i,1})
    sub_C = XY_Control(result{i,1}, result{i,2});
    sub_T = XY_Treat(result{i,1},   result{i,2});
    L = [L; mean(sub_C(:)), mean(sub_T(:))];
    i = i + 1;
end

fprintf('\n%-10s  %-12s  %-12s  %-16s  %-16s\n', ...
    'Cluster', 'X members', 'Y members', 'Control mean r', 'Treat mean r');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:n_clusters
    fprintf('%-10d  %-12d  %-12d  %-16.4f  %-16.4f\n', ...
        i, numel(result{i,1}), numel(result{i,2}), L(i,1), L(i,2));
end

%% Treat vs Control Bar Chart
figure;
bar(L, 'grouped');
set(gca, 'XTickLabel', arrayfun(@(i) sprintf('C%d',i), 1:n_clusters, 'UniformOutput', false));
legend({'Control', 'Treat'}, 'Box', 'off', 'Location', 'NorthEast');
ylabel('Mean Pearson r');
title('Per-cluster mean cross-correlation: Treat vs Control');
grid off;

%% t-statistic Distribution: Inside vs Outside Clusters
% Converts R_Control to a t-statistic and compares the distribution
% of values inside detected clusters vs the background (outside all clusters).
T = R_Control .* sqrt((127 - 2) ./ (1 - R.^2));

select_X = cat(2, result{1:n_clusters, 1});
select_Y = cat(2, result{1:n_clusters, 2});
omit_X   = setdiff(1:m, select_X);
omit_Y   = setdiff(1:n, select_Y);

data_out = T(omit_X, m + omit_Y);
data_out = data_out(:);

data_in = [];
for i = 1:2
    in = T(result{i,1}, m + result{i,2});
    data_in = [data_in; in(:)];
end

nbins = 50;
edges = linspace(min([data_out; data_in]), max([data_out; data_in]), nbins + 1);

figure; hold on;
histogram(data_out, edges, 'Normalization', 'pdf', ...
    'FaceColor', [0 0 0.7], 'FaceAlpha', 0.7, 'EdgeColor', 'k', 'LineWidth', 0.3);
histogram(data_in,  edges, 'Normalization', 'pdf', ...
    'FaceColor', [0 1 1],   'FaceAlpha', 0.5, 'EdgeColor', 'k', 'LineWidth', 0.3);
l1 = plot(NaN, NaN, '-', 'Color', [0 1 1],   'LineWidth', 2);
l2 = plot(NaN, NaN, '-', 'Color', [0 0 0.7], 'LineWidth', 2);
legend([l1 l2], {'In first two clusters', 'Outside clusters'}, ...
    'Interpreter', 'latex', 'Box', 'off', 'Location', 'NorthEast');
set(gca, 'FontSize', 20);
xlabel('t-statistic'); ylabel('Density');
title('t-statistic distribution: inside vs outside detected clusters');
grid off;

%% Detected Gene Lists
idx_X_all = cat(2, result{1:n_clusters, 1});
idx_Y_all = m + cat(2, result{1:n_clusters, 2});

name_X = name(idx_X_all);
name_Y = name(idx_Y_all);

for i = 1:n_clusters
    ix = result{i,1};
    iy = m + result{i,2};
    fprintf('\n=== Cluster %d  (%d X  x  %d Y genes) ===\n', ...
        i, numel(ix), numel(iy));
    fprintf('X genes: ');  fprintf('%s  ', name{ix});    fprintf('\n');
    fprintf('Y genes: ');  fprintf('%s  ', name{iy});    fprintf('\n');
end
