%% Dense Community Detection: CD vs UC vs nonIBD
%
% Identifies co-expression communities across three diagnostic groups:
%   CD      — Crohn's Disease
%   UC      — Ulcerative Colitis
%   nonIBD  — healthy controls
%
% Input: CleanCorNonzero.mat
%   R1, P1 — correlation matrices for CD group
%   R2, P2 — correlation matrices for UC group
%   R3, P3 — correlation matrices for nonIBD group
%   name   — feature names (m transcriptomics + n metabolomics)
%
% Pipeline:
% # Load pre-computed correlation matrices for all three groups
% # Build combined P and R: max significance, max absolute correlation
% # Construct bipartite weight matrix Wp, threshold at -log10(p) >= 2
% # Greedy peeling to detect dense bipartite communities
% # Refine Cluster 2 by re-peeling its X side
% # Sort and visualise results
% # Report per-cluster mean correlation across CD, UC, nonIBD

%% Setup
addpath('../../methods')
addpath('../../visualization')

set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultTextFontSize', 18);
set(0, 'DefaultAxesColor', 'white');
set(0, 'DefaultFigureColor', 'white');

%% Load Correlation Matrices
% R1/P1 = CD, R2/P2 = UC, R3/P3 = nonIBD
% m = 17026 transcriptomics features (X), n = 476 metabolomics features (Y)
load('CleanCorNonzero.mat');

m = 17026;
n = 476;

fprintf('Loaded CleanCorNonzero.mat\n');
fprintf('  X features (transcriptomics) : %d\n', m);
fprintf('  Y features (metabolomics)    : %d\n', n);

%% Build Combined Correlation and Significance Matrices
% P = element-wise max across all three groups — keeps strongest signal.
% R = element-wise max by absolute value across CD, UC, nonIBD.
P = max(max(P1, P2), P3);

[~, idx] = max(cat(3, abs(R1), abs(R2), abs(R3)), [], 3);
R = R1;
R(idx == 2) = R2(idx == 2);
R(idx == 3) = R3(idx == 3);

%% Construct Bipartite Weight Matrix
% Wp = P(X, Y) block, entries below -log10(0.01) = 2 zeroed out.
Wp = P(1:m, (m+1):(m+n));
Wp(Wp < 2) = 0;

figure;
imagesc(P((m+1):(m+n), 1:m)); colorbar; colormap(jet); caxis([0 10]);
title('Bipartite significance matrix W_p  (Y \times X)');
xlabel('X (transcriptomics)'); ylabel('Y (metabolomics)');

%% Initial Community Detection: Greedy Peeling
% Threshold 1.6 controls minimum subgraph density.
result = greedy_peeling_XY_all(Wp, 1.6);
result([3 4], :) = result([4 3], :);
result(4, :) = {[]};

%% Refine Cluster 2: Re-peel X Side
idx_X = cat(1, result{2,1});
X = R(idx_X, idx_X);
X = (X + X') / 2;
test = greedy_peeling_X_all(X, 1.9);
test = sort_result(X, test);
result{2,1} = idx_X(cat(2, test{:,1}));

%% Sort and Finalise Clusters
idx_X = 1:m;
idx_Y = (m+1):(m+n);
X  = R(idx_X, idx_X);
Y  = R(idx_Y, idx_Y);
XY = R(idx_X, idx_Y);
result = sort_result_XandY(X, Y, result);
result = result(:, [2 1]);

n_clusters = sum(~cellfun(@isempty, result(:,1)));
fprintf('Clusters detected: %d\n', n_clusters);

%% Heatmaps: Combined Significance and Correlation
plot_all_result(P(idx_Y, idx_X), result, false, true, [0 10]);

%% 3-Panel Overview: Combined R (CD-dominant)
plot3in1_v1(X, Y, XY, result, true, [200, 5], true);

%% 3-Panel Overview: nonIBD Group (R3)
X3  = R3(idx_X, idx_X);
Y3  = R3(idx_Y, idx_Y);
XY3 = R3(idx_X, idx_Y);
plot3in1_v1(X3, Y3, XY3, result, true, [200, 5], true);

%% Per-cluster Mean Correlation: CD vs UC vs nonIBD
% Rows = result clusters, columns = [CD mean, UC mean, nonIBD mean].
% Compares whether each community is specific to a diagnostic group.
XY1 = R1(idx_X, idx_Y);
XY2 = R2(idx_X, idx_Y);
XY3 = R3(idx_X, idx_Y);

fprintf('\n%-10s  %-12s  %-12s  %-12s  %-12s  %-12s\n', ...
    'Cluster', 'X members', 'Y members', 'CD mean r', 'UC mean r', 'nonIBD mean r');
fprintf('%s\n', repmat('-', 1, 74));

for i = 1:n_clusters
    rows = result{i, 1};   % metabolomics indices (Y, after column swap)
    cols = result{i, 2};   % transcriptomics indices (X)
    if ~isempty(rows) && ~isempty(cols)
        sub1 = XY1(cols, rows);
        sub2 = XY2(cols, rows);
        sub3 = XY3(cols, rows);
        fprintf('%-10d  %-12d  %-12d  %-12.4f  %-12.4f  %-12.4f\n', ...
            i, numel(rows), numel(cols), ...
            mean(sub1(:), 'omitnan'), mean(sub2(:), 'omitnan'), mean(sub3(:), 'omitnan'));
    end
end

%% Detected Gene / Metabolite Lists
result_named      = result;
result_named(:,1) = cellfun(@(x) x + m, result(:,1), 'UniformOutput', false);
NameList          = cellfun(@(ix) name(ix), result_named, 'UniformOutput', false);

for i = 1:n_clusters
    fprintf('\n=== Cluster %d ===\n', i);
    fprintf('Metabolomics (%d):\n',    numel(NameList{i,1}));
    fprintf('  %s\n', NameList{i,1}{:});
    fprintf('Transcriptomics (%d):\n', numel(NameList{i,2}));
    fprintf('  %s\n', NameList{i,2}{:});
end
