function plot_comparison(results_all, cfg, fig_dir)
% PLOT_COMPARISON  Generate comparison figures for FDI algorithm evaluation.
%   plot_comparison(results_all, cfg, fig_dir) creates two comparison figures:
%
%   Figure 1: Raw ISIS cross-check vs. SWGLT for hard fault conditions (1 & 2)
%     - Shows that naive ISIS residual detection has high false alarm rate
%     - Shows that SWGLT provides clean detection with low FAR
%
%   Figure 2: WGLT vs. SWGLT for soft/ramp fault (Condition 3)
%     - Shows WGLT (single-frame) has longer detection delay
%     - Shows SWGLT (sliding window) detects soft faults earlier
%
%   Inputs:
%     results_all - 1 x 3 struct array with fields per condition:
%                   .t, .FD_wglt, .FD_swglt, .FD_raw, .T_raw, .T_D,
%                   .T_adaptive, .fault_mask, .fault_intervals, .cond_title,
%                   .is_soft, .fault_cfg, .stats_wglt, .stats_swglt,
%                   .stats_raw
%     cfg         - configuration structure
%     fig_dir     - output directory for saved figures

    plot_raw_vs_algorithm(results_all, fig_dir);
    plot_soft_fault_comparison(results_all, cfg, fig_dir);
end


%% ========================================================================
%  Figure 1: Raw ISIS Cross-check vs. SWGLT Algorithm
% =========================================================================
function plot_raw_vs_algorithm(results_all, fig_dir)

    fig = figure('Name', 'Comparison: Raw ISIS vs SWGLT', ...
        'Position', [30 30 1400 900], 'Color', 'w');

    row_titles = {};
    for cond = 1:2
        R = results_all(cond);
        t = R.t;

        % Determine which INS has the fault
        if R.fault_cfg.sensor_idx <= 3
            ins_idx = 1;
            ins_label = 'INS1';
        else
            ins_idx = 2;
            ins_label = 'INS2';
        end

        row_titles{cond} = R.cond_title;

        % --- Left subplot: Raw ISIS detection ---
        ax1 = subplot(2, 2, (cond-1)*2 + 1);
        hold on;
        yl_raw = [0, max([max(R.FD_raw(:, ins_idx)), R.T_raw]) * 1.5];
        if yl_raw(2) <= yl_raw(1), yl_raw(2) = yl_raw(1) + 1; end

        draw_fault_patches(t, R.fault_intervals, yl_raw);

        plot(t, R.FD_raw(:, ins_idx), 'Color', [0.6 0.6 0.6], 'LineWidth', 0.3);
        hTh = yline(R.T_raw, 'r--', 'LineWidth', 1.5);

        xlim([t(1) t(end)]);
        ylim(yl_raw);
        ylabel('检测统计量');
        if cond == 2, xlabel('时间 (s)'); end
        title(sprintf('直接对比法 (%s vs ISIS)', ins_label), 'FontSize', 11);
        grid on;

        % Annotation box with FAR/FDR
        far_raw = R.stats_raw.FAR * 100;
        fdr_raw = R.stats_raw.FDR * 100;
        text_str = sprintf('FAR=%.1f%%\nFDR=%.1f%%', far_raw, fdr_raw);
        add_annotation_box(ax1, text_str, [0.6 0.6 0.6]);
        hold off;

        % --- Right subplot: SWGLT detection ---
        ax2 = subplot(2, 2, (cond-1)*2 + 2);
        hold on;
        yl_sw = [0, max([max(R.FD_swglt), R.T_adaptive]) * 1.3];
        if yl_sw(2) <= yl_sw(1), yl_sw(2) = yl_sw(1) + 1; end

        draw_fault_patches(t, R.fault_intervals, yl_sw);

        plot(t, R.FD_swglt, 'Color', [0.1 0.7 0.3], 'LineWidth', 0.5);
        yline(R.T_adaptive, 'r--', 'LineWidth', 1.5);

        xlim([t(1) t(end)]);
        ylim(yl_sw);
        ylabel('检测统计量');
        if cond == 2, xlabel('时间 (s)'); end
        title('SWGLT算法 (滑动窗口加权奇偶向量法)', 'FontSize', 11);
        grid on;

        % Annotation box with FAR/FDR
        far_sw = R.stats_swglt.FAR * 100;
        fdr_sw = R.stats_swglt.FDR * 100;
        text_str = sprintf('FAR=%.1f%%\nFDR=%.1f%%', far_sw, fdr_sw);
        add_annotation_box(ax2, text_str, [0.1 0.7 0.3]);
        hold off;
    end

    sgtitle({'未处理ISIS信号 vs 优化FDI算法检测效果对比', ''}, ...
        'FontSize', 14, 'FontWeight', 'bold');

    saveas(fig, fullfile(fig_dir, 'comparison_raw_vs_algorithm.fig'));
    saveas(fig, fullfile(fig_dir, 'comparison_raw_vs_algorithm.png'));
    fprintf('Saved: comparison_raw_vs_algorithm.fig/.png\n');
end


%% ========================================================================
%  Figure 2: WGLT vs. SWGLT for Soft Fault Detection
% =========================================================================
function plot_soft_fault_comparison(results_all, cfg, fig_dir)

    R = results_all(3);  % Condition 3 = soft fault
    if ~R.is_soft
        warning('Condition 3 is not marked as soft fault. Skipping comparison plot.');
        return;
    end

    t = R.t;
    fault_start = R.fault_intervals(1);
    fault_end   = R.fault_intervals(2);

    % Zoom window: 20s before fault start to 20s after fault end
    t_zoom_start = max(fault_start - 20, t(1));
    t_zoom_end   = min(fault_end + 20, t(end));
    zoom_mask = (t >= t_zoom_start) & (t <= t_zoom_end);
    t_zoom = t(zoom_mask);

    fig = figure('Name', 'Comparison: Soft Fault - WGLT vs SWGLT', ...
        'Position', [60 60 1200 700], 'Color', 'w');

    % --- Top subplot: WGLT (original algorithm) ---
    ax1 = subplot(2, 1, 1);
    hold on;
    FD_wglt_zoom = R.FD_wglt(zoom_mask);
    yl_w = [0, max([max(FD_wglt_zoom), R.T_D]) * 1.5];
    if yl_w(2) <= yl_w(1), yl_w(2) = yl_w(1) + 1; end

    % Fault region shading
    patch([fault_start fault_end fault_end fault_start], ...
        [yl_w(1) yl_w(1) yl_w(2) yl_w(2)], ...
        [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    plot(t_zoom, FD_wglt_zoom, 'Color', [0.8 0.5 0.1], 'LineWidth', 0.8);
    yline(R.T_D, 'r--', 'LineWidth', 1.5);

    % Mark fault start
    xline(fault_start, 'k:', 'LineWidth', 1.0, 'Label', '故障起始');

    % Mark WGLT detection time
    delay_wglt = R.stats_wglt.delay;
    if ~isnan(delay_wglt) && ~isinf(delay_wglt)
        t_det_wglt = fault_start + delay_wglt;
        xline(t_det_wglt, 'b-', 'LineWidth', 1.5, 'Label', ...
            sprintf('检测时刻 (延迟%.1fs)', delay_wglt));
    end

    xlim([t_zoom_start t_zoom_end]);
    ylim(yl_w);
    ylabel('检测统计量');
    title(sprintf('WGLT (单帧检测, 原始算法)  —  FDR=%.1f%%, FAR=%.1f%%', ...
        R.stats_wglt.FDR*100, R.stats_wglt.FAR*100), 'FontSize', 12);
    grid on;
    hold off;

    % --- Bottom subplot: SWGLT (improved algorithm) ---
    ax2 = subplot(2, 1, 2);
    hold on;
    FD_swglt_zoom = R.FD_swglt(zoom_mask);
    yl_s = [0, max([max(FD_swglt_zoom), R.T_adaptive]) * 1.5];
    if yl_s(2) <= yl_s(1), yl_s(2) = yl_s(1) + 1; end

    % Fault region shading
    patch([fault_start fault_end fault_end fault_start], ...
        [yl_s(1) yl_s(1) yl_s(2) yl_s(2)], ...
        [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    plot(t_zoom, FD_swglt_zoom, 'Color', [0.1 0.7 0.3], 'LineWidth', 0.8);
    yline(R.T_adaptive, 'r--', 'LineWidth', 1.5);

    % Mark fault start
    xline(fault_start, 'k:', 'LineWidth', 1.0, 'Label', '故障起始');

    % Mark SWGLT detection time
    delay_swglt = R.stats_swglt.delay;
    if ~isnan(delay_swglt) && ~isinf(delay_swglt)
        t_det_swglt = fault_start + delay_swglt;
        xline(t_det_swglt, 'b-', 'LineWidth', 1.5, 'Label', ...
            sprintf('检测时刻 (延迟%.1fs)', delay_swglt));
    end

    xlim([t_zoom_start t_zoom_end]);
    ylim(yl_s);
    ylabel('检测统计量');
    xlabel('时间 (s)');
    title(sprintf('SWGLT (滑动窗口检测, 改进算法)  —  FDR=%.1f%%, FAR=%.1f%%', ...
        R.stats_swglt.FDR*100, R.stats_swglt.FAR*100), 'FontSize', 12);
    grid on;
    hold off;

    sgtitle({'缓变故障检测对比: 原始算法 vs 改进算法', ...
             sprintf('故障区间 [%.0f, %.0f]s', fault_start, fault_end)}, ...
        'FontSize', 14, 'FontWeight', 'bold');

    % Add delay comparison textbox
    if ~isnan(delay_wglt) && ~isinf(delay_wglt) && ~isnan(delay_swglt) && ~isinf(delay_swglt)
        improvement = delay_wglt - delay_swglt;
        if improvement > 0
            msg = sprintf('检测延迟改善: %.1fs -> %.1fs (提前 %.1fs)', ...
                delay_wglt, delay_swglt, improvement);
        else
            msg = sprintf('WGLT延迟: %.1fs, SWGLT延迟: %.1fs', delay_wglt, delay_swglt);
        end
        annotation('textbox', [0.15 0.01 0.7 0.04], 'String', msg, ...
            'HorizontalAlignment', 'center', 'FontSize', 11, ...
            'FontWeight', 'bold', 'EdgeColor', 'none', 'Color', [0 0.4 0]);
    end

    saveas(fig, fullfile(fig_dir, 'comparison_soft_fault.fig'));
    saveas(fig, fullfile(fig_dir, 'comparison_soft_fault.png'));
    fprintf('Saved: comparison_soft_fault.fig/.png\n');
end


%% ========================================================================
%  Helper functions
% =========================================================================
function draw_fault_patches(t, fault_intervals, yl)
    for fi = 1:size(fault_intervals, 1)
        patch([fault_intervals(fi,1) fault_intervals(fi,2) ...
               fault_intervals(fi,2) fault_intervals(fi,1)], ...
              [yl(1) yl(1) yl(2) yl(2)], ...
              [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
              'HandleVisibility', 'off');
    end
end


function add_annotation_box(ax, text_str, color)
    % Add a text annotation box in the top-left corner of the axes
    pos = get(ax, 'Position');
    annotation('textbox', [pos(1)+0.01, pos(2)+pos(4)-0.08, 0.12, 0.07], ...
        'String', text_str, 'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', color, 'BackgroundColor', [1 1 1 0.85], ...
        'EdgeColor', color, 'LineWidth', 1.2, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
