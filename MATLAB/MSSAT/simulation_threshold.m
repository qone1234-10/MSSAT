%% =========================================================
%% MSSAT sequential test statistic histograms + thresholds
%% using explicit fields from jive_multi_2:
%%   standardized statistic = out.test_path.Z
%%   threshold              = out.test_path.crit
%% =========================================================

var_keep = {'result_real', 'finalTable'};
vars = setdiff(who, var_keep);
clear(vars{:});
clc

addpath('C:\Users\qone\Desktop\Matlab\MyJIVE');
addpath('C:\Users\qone\Desktop\Matlab\AJIVE_Project\AJIVECode');
addpath('C:\Users\qone\Desktop\Matlab\Data-Integration-Via-Analysis-of-Subspaces\DJIVECode');

%% ---------------------------------------------------------
%% simulation setting
%% ---------------------------------------------------------
n = 500;
p1 = 400;
p2 = 300;
r1 = 6;
r2 = 9;
r  = 3;
d = 1;
ang = pi/8;

% d1 = d/2 * 3 * sqrt(n) * ((r1:-1:1)+1) ;
% d2 = d/2 * 2 * sqrt(n) * ((r2:-1:1)+1) ;

d1 = d/2 * 3 * sqrt(n) * [6,5,4,2,3,7] ;
d2 = d/2 * 2 * sqrt(n) * [7,6,5,2,3,4,8,9,10];

% d1 = d/4 * 4 * sqrt(n) * ((r1:-1:1)+1) ;
% d2 = d/4 * 3 * sqrt(n) * ((r2:-1:1)+1) ;

s1 = 1;
s2 = 1;
rep = 100;

v0 = random_orthonormal(n, r);
v_temp = random_orthonormal_orthogonal(r1 + r2 - 2 * r, v0);
v1 = v_temp(:, 1:(r1 - r));
v2 = [ ...
    v_temp(:, (r1-r+1):(2*r1 - 2*r)) * sin(ang) + v1 * cos(ang), ...
    v_temp(:, ((r1 - r)*2+1):end) ...
];

u1 = random_orthonormal(p1, r1);
u2 = random_orthonormal(p2, r2);

%% ---------------------------------------------------------
%% containers
%% ---------------------------------------------------------
ntest = 4;

z_alpha = nan(rep, ntest);
z_alt5  = nan(rep, ntest);
z_alt10 = nan(rep, ntest);

thr_1 = nan(rep, ntest);
thr_2 = nan(rep, ntest);
thr_3 = nan(rep, ntest);
thr_4 = nan(rep, ntest);
thr_5 = nan(rep, ntest);

%% ---------------------------------------------------------
%% simulation loop
%% ---------------------------------------------------------
rng(1);

for rep_idx = 1:rep

    rot = random_orthonormal(r, r);

    X1 = u1 * diag(d1) * [v0, v1]'     + s1 * randn(p1, n);
    X2 = u2 * diag(d2) * [v0*rot, v2]' + s2 * randn(p2, n);
    
    [s1_hat2, r1_hat] = BEMA_combined(X1);
    [s2_hat2, r2_hat] = BEMA_combined(X2);

    X_list = {X1, X2};

    %% ---------------- method 1: alpha = 0.05 ----------------
    out_alpha = jive_multi_3( ...
        X_list, ...
        'center',       false, ...
        'alpha',        0.05, ...
        'rank_method',  'given', ...
        'sigma_method', 'given', ...
        'r_init',       [r1_hat, r2_hat], ...
        'sigma_init',   [s1_hat2, s2_hat2], ...
        'test_mode',    'alpha_seq', ...
        'mode',         'test', ...
        'alt_angle_deg', 5, ...
        'verbose',      false ...
    );

    %% ---------------- method 4: altmean 5 degree ----------------
    out_alt5 = jive_multi_3( ...
        X_list, ...
        'center',       false, ...
        'alpha',        0.05, ...
        'rank_method',  'given', ...
        'sigma_method', 'given', ...
        'r_init',       [r1_hat, r2_hat], ...
        'sigma_init',   [s1_hat2, s2_hat2], ...
        'test_mode',    'altmean', ...
        'mode',         'test', ...
        'alt_angle_deg', 5, ...
        'verbose',      false ...
    );

    %% ---------------- method 5: altmean 10 degree ----------------
    out_alt10 = jive_multi_3( ...
        X_list, ...
        'center',       false, ...
        'alpha',        0.05, ...
        'rank_method',  'given', ...
        'sigma_method', 'given', ...
        'r_init',       [r1_hat, r2_hat], ...
        'sigma_init',   [s1_hat2, s2_hat2], ...
        'test_mode',    'altmean', ...
        'mode',         'test', ...
        'alt_angle_deg', 10, ...
        'verbose',      false ...
    );

    %% standardized statistics from explicit fields
    ztmp = out_alpha.test_path.Z(:)';
    z_alpha(rep_idx, 1:min(ntest, numel(ztmp))) = ztmp(1:min(ntest, numel(ztmp)));

    ztmp = out_alt5.test_path.Z(:)';
    z_alt5(rep_idx, 1:min(ntest, numel(ztmp))) = ztmp(1:min(ntest, numel(ztmp)));

    ztmp = out_alt10.test_path.Z(:)';
    z_alt10(rep_idx, 1:min(ntest, numel(ztmp))) = ztmp(1:min(ntest, numel(ztmp)));

    %% thresholds: methods 1,2,3
    for s = 1:ntest
        thr_1(rep_idx, s) = norminv(0.05);          % one-sided lower-tail cutoff on Z
        thr_2(rep_idx, s) = norminv(0.05 / (2^s));  % sequential Bonferroni-style
        thr_3(rep_idx, s) = -6;                     % 6 sigma lower-tail cutoff
    end

    %% thresholds: methods 4,5 from explicit fields
    mu0_tmp   = out_alt5.test_path.mu0(:)';
    muAlt_tmp = out_alt5.test_path.mu_alt(:)';
    V0_tmp    = out_alt5.test_path.V0(:)';

    m = min(ntest, min([numel(mu0_tmp), numel(muAlt_tmp), numel(V0_tmp)]));
    thr_4(rep_idx, 1:m) = sqrt(n ./ V0_tmp(1:m)) .* (muAlt_tmp(1:m) - mu0_tmp(1:m));

    mu0_tmp   = out_alt10.test_path.mu0(:)';
    muAlt_tmp = out_alt10.test_path.mu_alt(:)';
    V0_tmp    = out_alt10.test_path.V0(:)';

    m = min(ntest, min([numel(mu0_tmp), numel(muAlt_tmp), numel(V0_tmp)]));
    thr_5(rep_idx, 1:m) = sqrt(n ./ V0_tmp(1:m)) .* (muAlt_tmp(1:m) - mu0_tmp(1:m));

    if mod(rep_idx, 10) == 0
        fprintf('rep %d / %d done\n', rep_idx, rep);
    end
end

%% ---------------------------------------------------------
%% mean thresholds for plotting
%% ---------------------------------------------------------
thr1_plot = mean(thr_1, 1, 'omitnan');
thr2_plot = mean(thr_2, 1, 'omitnan');
thr3_plot = mean(thr_3, 1, 'omitnan');
thr4_plot = mean(thr_4, 1, 'omitnan');
thr5_plot = mean(thr_5, 1, 'omitnan');

%% ---------------------------------------------------------
%% 3 plots = sequential tests 1,2,3
% %% ---------------------------------------------------------
% figure('Color','w','Position',[100 100 1600 420]);
% tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');
% 
% for s = 1:ntest
%     nexttile;
%     hold on;
% 
%     % histograms of standardized statistics
%     histogram(z_alpha(:,s), ...
%         'Normalization','pdf', ...
%         'FaceColor',[0.2 0.4 0.8], ...
%         'FaceAlpha',0.20, ...
%         'EdgeColor','none');
% 
%     histogram(z_alt5(:,s), ...
%         'Normalization','pdf', ...
%         'FaceColor',[0.85 0.25 0.25], ...
%         'FaceAlpha',0.20, ...
%         'EdgeColor','none');
% 
%     histogram(z_alt10(:,s), ...
%         'Normalization','pdf', ...
%         'FaceColor',[0.85 0.55 0.25], ...
%         'FaceAlpha',0.20, ...
%         'EdgeColor','none');
% 
%     % threshold lines
%     xline(thr1_plot(s), '-', ...
%         'Color',[0 0.2 0.8], 'LineWidth',2.0);
% 
%     xline(thr2_plot(s), '--', ...
%         'Color',[0 0.2 0.8], 'LineWidth',2.0);
% 
%     xline(thr3_plot(s), '-.', ...
%         'Color',[0 0.2 0.8], 'LineWidth',2.0);
% 
%     xline(thr4_plot(s), '-', ...
%         'Color',[0.8 0 0], 'LineWidth',2.0);
% 
%     xline(thr5_plot(s), '--', ...
%         'Color',[0.8 0 0], 'LineWidth',2.0);
% 
%     xlabel('Standardized test statistic');
%     ylabel('Empirical density');
%     title(sprintf('Sequential test %d', s));
%     grid on;
%     box on;
% 
%     hold off;
% end
% 
% sgtitle('Empirical histograms of sequential standardized test statistics with thresholds');
% 
% legend( ...
%     {'alpha-mode stat', 'altmean 5° stat', 'altmean 10° stat', ...
%      '0.05', '0.05 / 2^s', '6\sigma', 'alt 5°', 'alt 10°'}, ...
%     'Location', 'bestoutside');

%% ---------------------------------------------------------
%% summary table
%% ---------------------------------------------------------
threshold_table = table( ...
    (1:ntest)', ...
    thr1_plot(:), ...
    thr2_plot(:), ...
    thr3_plot(:), ...
    thr4_plot(:), ...
    thr5_plot(:), ...
    'VariableNames', {'seq_test','thr_0p05','thr_0p05_over_2pow_s','thr_6sigma','thr_alt5','thr_alt10'} ...
);

disp('===== Threshold summary =====');
disp(threshold_table);

%% ---------------------------------------------------------
%% save in workspace
%% ---------------------------------------------------------
seq_test_hist_result = struct();
seq_test_hist_result.z_alpha = z_alpha;
seq_test_hist_result.z_alt5  = z_alt5;
seq_test_hist_result.z_alt10 = z_alt10;

seq_test_hist_result.thr_1 = thr_1;
seq_test_hist_result.thr_2 = thr_2;
seq_test_hist_result.thr_3 = thr_3;
seq_test_hist_result.thr_4 = thr_4;
seq_test_hist_result.thr_5 = thr_5;

seq_test_hist_result.thr_plot = threshold_table;

disp('Saved result struct: seq_test_hist_result');


%% ---------------------------------------------------------
%% Broken x-axis histogram figure: 4 sequential tests, 2 by 2
%% methods 1,2,3: histogram of standardized test statistic + threshold lines
%% methods 4,5: histogram of standardized thresholds
%% ---------------------------------------------------------
% 
% ntest = 4;   % plot sequential tests 1,2,3,4
% 
% figure('Color','w','Position',[80 80 1450 900]);
% 
% % 각 sequential test별 x-range 수동 지정
% % 필요하면 결과 보고 숫자만 조정하면 됨
% xlims_left  = {[-110 -70], [-95 -55], [-75 -40], [-65 -35]};
% xlims_right = {[-25 5],    [-25 5],   [-20 5],   [-20 5]};
% 
% % colors
% col_stat  = [0.45 0.45 0.45];   % 회색
% col_alt5  = [1.00 0.50 0.50];   % 매우 연한 빨간색
% col_alt10 = [1.00 0.20 0.20];   % 연한 빨간색
% col_line  = [0.00 0.20 0.80];   % 파란 threshold line
% 
% % title
% annotation('textbox', [0.22 0.955 0.56 0.04], ...
%     'String', 'Empirical histograms of sequential standardized statistics / thresholds', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center', ...
%     'FontWeight', 'bold', ...
%     'FontSize', 14);
% 
% for s = 1:ntest
% 
%     outer = subplot(2,2,s);
%     outer_pos = get(outer, 'Position');
%     delete(outer);
% 
%     left_ratio  = 0.22;
%     gap_ratio   = 0.04;
%     right_ratio = 1 - left_ratio - gap_ratio;
% 
%     % 오른쪽 legend 공간 확보: 오른쪽 column만 약간 줄임
%     if mod(s, 2) == 0
%         outer_pos(3) = outer_pos(3) * 0.82;
%     else
%         outer_pos(3) = outer_pos(3) * 0.92;
%     end
% 
%     axL = axes('Position', [outer_pos(1), ...
%                             outer_pos(2), ...
%                             outer_pos(3)*left_ratio, ...
%                             outer_pos(4)]);
% 
%     axR = axes('Position', [outer_pos(1)+outer_pos(3)*(left_ratio+gap_ratio), ...
%                             outer_pos(2), ...
%                             outer_pos(3)*right_ratio, ...
%                             outer_pos(4)]);
% 
%     %% ---------------- left axis ----------------
%     axes(axL); hold on;
% 
%     h1 = histogram(z_alpha(:,s), ...
%         'Normalization', 'pdf', ...
%         'FaceColor', col_stat, ...
%         'FaceAlpha', 0.40, ...
%         'EdgeColor', 'none');
% 
%     h2 = histogram(thr_4(:,s), ...
%         'Normalization', 'pdf', ...
%         'FaceColor', col_alt5, ...
%         'FaceAlpha', 0.60, ...
%         'EdgeColor', 'none');
% 
%     h3 = histogram(thr_5(:,s), ...
%         'Normalization', 'pdf', ...
%         'FaceColor', col_alt10, ...
%         'FaceAlpha', 0.55, ...
%         'EdgeColor', 'none');
% 
%     l1 = xline(thr1_plot(s), '-', ...
%         'Color', col_line, 'LineWidth', 2.0);
% 
%     l2 = xline(thr2_plot(s), '--', ...
%         'Color', col_line, 'LineWidth', 2.0);
% 
%     l3 = xline(thr3_plot(s), '-.', ...
%         'Color', col_line, 'LineWidth', 2.0);
% 
%     xlim(xlims_left{s});
%     grid on;
%     box on;
%     ylabel('Empirical density');
%     title(sprintf('Sequential test %d', s));
% 
%     %% ---------------- right axis ----------------
%     axes(axR); hold on;
% 
%     histogram(z_alpha(:,s), ...
%         'Normalization', 'pdf', ...
%         'FaceColor', col_stat, ...
%         'FaceAlpha', 0.40, ...
%         'EdgeColor', 'none');
% 
%     histogram(thr_4(:,s), ...
%         'Normalization', 'pdf', ...
%         'FaceColor', col_alt5, ...
%         'FaceAlpha', 0.60, ...
%         'EdgeColor', 'none');
% 
%     histogram(thr_5(:,s), ...
%         'Normalization', 'pdf', ...
%         'FaceColor', col_alt10, ...
%         'FaceAlpha', 0.55, ...
%         'EdgeColor', 'none');
% 
%     xline(thr1_plot(s), '-', ...
%         'Color', col_line, 'LineWidth', 2.0);
% 
%     xline(thr2_plot(s), '--', ...
%         'Color', col_line, 'LineWidth', 2.0);
% 
%     xline(thr3_plot(s), '-.', ...
%         'Color', col_line, 'LineWidth', 2.0);
% 
%     xlim(xlims_right{s});
%     grid on;
%     box on;
%     xlabel('Standardized scale');
%     set(axR, 'YTickLabel', []);
% 
%     % y-limit 통일
%     ylL = ylim(axL);
%     ylR = ylim(axR);
%     ymax = max([ylL(2), ylR(2)]);
%     ylim(axL, [0 ymax]);
%     ylim(axR, [0 ymax]);
% 
%     % style
%     axL.Box = 'off';
%     axR.Box = 'off';
%     axL.YAxisLocation = 'left';
%     axR.YAxisLocation = 'right';
%     axR.YTick = [];
%     axL.XColor = [0 0 0];
%     axR.XColor = [0 0 0];
%     axL.YColor = [0 0 0];
%     axR.YColor = [0 0 0];
% 
%     % break mark //
%     annotation('line', ...
%         [axL.Position(1)+axL.Position(3)-0.004, axL.Position(1)+axL.Position(3)+0.004], ...
%         [axL.Position(2)+0.02, axL.Position(2)+0.04], ...
%         'LineWidth', 1.5, 'Color', 'k');
%     annotation('line', ...
%         [axL.Position(1)+axL.Position(3)+0.006, axL.Position(1)+axL.Position(3)+0.014], ...
%         [axL.Position(2)+0.02, axL.Position(2)+0.04], ...
%         'LineWidth', 1.5, 'Color', 'k');
% 
%     annotation('line', ...
%         [axL.Position(1)+axL.Position(3)-0.004, axL.Position(1)+axL.Position(3)+0.004], ...
%         [axL.Position(2)+axL.Position(4)-0.04, axL.Position(2)+axL.Position(4)-0.02], ...
%         'LineWidth', 1.5, 'Color', 'k');
%     annotation('line', ...
%         [axL.Position(1)+axL.Position(3)+0.006, axL.Position(1)+axL.Position(3)+0.014], ...
%         [axL.Position(2)+axL.Position(4)-0.04, axL.Position(2)+axL.Position(4)-0.02], ...
%         'LineWidth', 1.5, 'Color', 'k');
% end
% 
% % legend
% lgd = legend([h1, h2, h3, l1, l2, l3], ...
%     {'standardized test stat', 'alt 5° threshold', 'alt 10° threshold', ...
%      '0.05', '0.05 / 2^s', '6\sigma'}, ...
%     'Location', 'eastoutside', ...
%     'FontSize', 10);
% 
% lgd.Position(1) = 0.89;
% lgd.Position(2) = 0.43;

%% =========================================================
%% Histogram figure without broken axis
%% =========================================================

ntest = 4;

figure('Color','w','Position',[80 80 1150 900]);

% x/y range manually 지정
xlims = {
    [-12 3.5], ...
    [-12 3.5], ...
    [-12 3.5], ...
    [-25 3.5]
};

% ylims = {
%     [0 2.70], ...
%     [0 2.70], ...
%     [0 2.70], ...
%     [0 2.70]
% };

ylims = {
    [0 1.00], ...
    [0 1.00], ...
    [0 1.00], ...
    [0 1.00]
};

% colors
col_stat  = [0.15 0.15 0.15];
col_alt5  = [1.00 0.60 0.60];
col_alt10 = [1.00 0.20 0.20];
col_line  = [0.00 0.40 0.90];

annotation('textbox', [0.22 0.955 0.56 0.04], ...
    'String', 'Empirical histograms of sequential standardized statistics / thresholds', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', ...
    'FontSize', 14);

for s = 1:ntest

    subplot(2,2,s);
    
    ax = gca;

    pos = ax.Position;

    % 오른쪽 여백 확보 (legend 공간)
    pos(3) = pos(3) * 0.75;

    ax.Position = pos;

    hold on;

    h1 = histogram(z_alpha(:,s), ...
        'Normalization', 'pdf', ...
        'FaceColor', col_stat, ...
        'FaceAlpha', 0.40, ...
        'EdgeColor', 'none');

    h2 = histogram(thr_4(:,s), ...
        'Normalization', 'pdf', ...
        'FaceColor', col_alt5, ...
        'FaceAlpha', 0.60, ...
        'EdgeColor', 'none');

    h3 = histogram(thr_5(:,s), ...
        'Normalization', 'pdf', ...
        'FaceColor', col_alt10, ...
        'FaceAlpha', 0.55, ...
        'EdgeColor', 'none');

    l1 = xline(thr1_plot(s), '-', ...
        'Color', col_line, 'LineWidth', 2.0);

    l2 = xline(thr2_plot(s), '--', ...
        'Color', col_line, 'LineWidth', 2.0);

    l3 = xline(thr3_plot(s), '-.', ...
        'Color', col_line, 'LineWidth', 2.0);

    % manually 지정
    xlim(xlims{s});
    ylim(ylims{s});

    grid on;
    box on;

    xlabel('Standardized scale');
    ylabel('Empirical density');

    title(sprintf('Sequential test %d', s));
    %% =====================================================
    %% y-axis break mark (fixed)
    %% =====================================================

    ax = gca;

    xl = xlim;

    % break 위치
    yb = 0.95;

    xr = xl(2) - xl(1);

    % ===============================
    % 1. 가로선 길이 (절반 수준)
    % ===============================
    dx = 0.025 * xr;   % ← 기존보다 확 줄임

    % ===============================
    % 2. y축 "끊김"용 흰 박스 (중요)
    % ===============================
    % rectangle( ...
    %     'Position', ...
    %     [xl(1)-0.25*xr, yb-0.06, 0.50*xr, 0.12], ...  % ← 충분히 크게 덮음
    %     'FaceColor', 'w', ...
    %     'EdgeColor', 'none', ...
    %     'Clipping', 'off');

    % ===============================
    % 3. 위/아래 가로선 (얇고 짧게)
    % ===============================
    lw = 1.2;   % ← 기존 2.4 → 0.6배 수준

    plot([xl(1)-dx, xl(1)+dx], ...
        [yb+0.02 yb+0.02], ...
        'k', 'LineWidth', lw, 'Clipping', 'off');

    plot([xl(1)-dx, xl(1)+dx], ...
        [yb-0.02 yb-0.02], ...
        'k', 'LineWidth', lw, 'Clipping', 'off');

    % ===============================
    % 4. tick (위쪽 제거 느낌)
    % ===============================
    yticks([0 0.2 0.4 0.6 0.8]);

    ax.Layer = 'top';
end

legend
lgd = legend([h1, h2, h3, l1, l2, l3], ...
    {'standardized test stat', ...
     'alt 5° threshold', ...
     'alt 10° threshold', ...
     '0.05', ...
     '0.05 / 2^s', ...
     '6\sigma'}, ...
    'Location', 'eastoutside', ...
    'FontSize', 10);

lgd.Position(1) = 0.83;
lgd.Position(2) = 0.46;

%% =====================================================
%% save figure (clean + publication quality)
%% =====================================================

out_name = 'seq_global_test_2block_2by2_broken_axis';

% PNG (high resolution)
exportgraphics(gcf, [out_name '.png'], ...
    'Resolution', 300);

% PDF (vector, 논문용 추천)
exportgraphics(gcf, [out_name '.pdf'], ...
    'ContentType', 'vector');

% FIG (MATLAB 재열기용)
savefig(gcf, [out_name '.fig']);

disp(['Saved figure: ' out_name]);

%% =====================================================
%% =====================================================
%% =====================================================
%% =====================================================
%% =====================================================
%% =====================================================


ntest = 4;

figure('Color','w','Position',[80 80 1150 900]);

% =========================
% tiled layout (핵심 변경)
% =========================
t = tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% 전체를 아래/왼쪽으로 밀기 + 오른쪽 공간 확보
t.OuterPosition = [0.03 0.03 0.85 0.92];

% x/y range manually 지정
xlims = {
    [-12 3.5], ...
    [-12 3.5], ...
    [-12 3.5], ...
    [-25 3.5]
};

ylims = {
    [0 1.00], ...
    [0 1.00], ...
    [0 1.00], ...
    [0 1.00]
};

% colors
col_stat  = [0.15 0.15 0.15];
col_alt5  = [1.00 0.60 0.60];
col_alt10 = [1.00 0.20 0.20];
col_line  = [0.00 0.40 0.90];

% title
annotation('textbox', [0.22 0.940 0.56 0.04], ...
    'String', 'Empirical histograms of sequential standardized statistics / thresholds', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', ...
    'FontSize', 14);

for s = 1:ntest

    ax = nexttile;
    ax.Position(2) = ax.Position(2) - 0.01;   % ↓ 살짝 아래로
    hold(ax,'on');

    h1 = histogram(ax, z_alpha(:,s), ...
        'Normalization', 'pdf', ...
        'FaceColor', col_stat, ...
        'FaceAlpha', 0.40, ...
        'EdgeColor', 'none');

    h2 = histogram(ax, thr_4(:,s), ...
        'Normalization', 'pdf', ...
        'FaceColor', col_alt5, ...
        'FaceAlpha', 0.60, ...
        'EdgeColor', 'none');

    h3 = histogram(ax, thr_5(:,s), ...
        'Normalization', 'pdf', ...
        'FaceColor', col_alt10, ...
        'FaceAlpha', 0.55, ...
        'EdgeColor', 'none');

    l1 = xline(ax, thr1_plot(s), '-', ...
        'Color', col_line, 'LineWidth', 2.0);

    l2 = xline(ax, thr2_plot(s), '--', ...
        'Color', col_line, 'LineWidth', 2.0);

    l3 = xline(ax, thr3_plot(s), '-.', ...
        'Color', col_line, 'LineWidth', 2.0);

    % limits
    xlim(ax, xlims{s});
    ylim(ax, ylims{s});

    grid(ax,'on');
    box(ax,'on');

    xlabel(ax,'Standardized scale');
    ylabel(ax,'Empirical density');
    title(ax, sprintf('Sequential test %d', s));

    %% =========================
    %% y-axis break mark
    %% =========================

    xl = xlim(ax);
    yb = 0.95;

    xr = xl(2) - xl(1);

    dx = 0.025 * xr;
    lw = 1.2;

    plot(ax, [xl(1)-dx, xl(1)+dx], [yb+0.02 yb+0.02], ...
        'k','LineWidth',lw,'Clipping','off');

    plot(ax, [xl(1)-dx, xl(1)+dx], [yb-0.02 yb-0.02], ...
        'k','LineWidth',lw,'Clipping','off');

    yticks(ax,[0 0.2 0.4 0.6 0.8]);

    ax.Layer = 'top';
end

lgd = legend([h1, h2, h3, l1, l2, l3], ...
    {'std stat', ...
     'alt 5°', ...
     'alt 10°', ...
     '0.05', ...
     '0.05/2^s', ...
     '6\sigma'}, ...
    'Location','eastoutside', ...
    'FontSize',8);   % ← 핵심 (줄임)

lgd.Box = 'off';

lgd.Position(1) = 0.95;
lgd.Position(2) = 0.46;

% =========================
% save
% =========================

out_name = 'seq_global_test_2block_2by2_tiledlayout';

exportgraphics(gcf, [out_name '.png'], 'Resolution', 300);
exportgraphics(gcf, [out_name '.pdf'], 'ContentType','vector');
savefig(gcf, [out_name '.fig']);

disp(['Saved figure: ' out_name]);