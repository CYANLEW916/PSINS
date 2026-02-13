function plot_comprehensive_comparison(all_results, cfg)
% PLOT_COMPREHENSIVE_COMPARISON  Cross-condition algorithm comparison.
%   Generates six comprehensive figures comparing GLT, WGLT, and SWGLT
%   across all fault conditions:
%     Figure 1: Performance metrics grouped bar charts (FDR, FAR, Accuracy)
%     Figure 2: Performance heatmap
%     Figure 3: Radar chart for multi-dimensional algorithm comparison
%     Figure 4: ROC curves (FDR vs FAR at varying thresholds)
%     Figure 5: Detection SNR and noise suppression analysis
%     Figure 6: Comprehensive summary dashboard
%
%   Input:
%     all_results - 1x3 struct array, each with fields:
%       .FD_glt, .FD_wglt, .FD_swglt  - detection functions (Nx1)
%       .FI_glt, .FI_wglt, .FI_swglt  - isolation functions (Nx9)
%       .T_D, .T_adaptive             - thresholds
%       .fault_mask                    - Nx1 logical
%       .stats_glt, .stats_wglt, .stats_swglt - performance structs
%       .t                            - Nx1 time vector
%       .cond_title                   - string
%     cfg - configuration structure from config.m
%
% See also  main_simulation, plot_results, evaluate_performance.

    methods = {'GLT', 'WGLT', 'SWGLT'};
    cond_labels = {'Cond 1: Gyro Hard', 'Cond 2: Accel Hard', 'Cond 3: Gyro Soft'};
    cond_labels_short = {'C1', 'C2', 'C3'};

    color_glt  = [0.55 0.55 0.55];
    color_wglt = [0.20 0.47 0.87];
    color_swglt= [0.85 0.22 0.16];
    colors_m = {color_glt, color_wglt, color_swglt};

    n_cond = length(all_results);

    % Extract all stats into matrices for easy access
    FDR_mat = zeros(n_cond, 3);  % conditions x methods
    FAR_mat = zeros(n_cond, 3);
    ACC_mat = zeros(n_cond, 3);
    DEL_mat = nan(n_cond, 3);

    for c = 1:n_cond
        stats = {all_results(c).stats_glt, all_results(c).stats_wglt, all_results(c).stats_swglt};
        for m = 1:3
            FDR_mat(c, m) = stats{m}.FDR * 100;
            FAR_mat(c, m) = stats{m}.FAR * 100;
            ACC_mat(c, m) = stats{m}.Accuracy * 100;
            if ~isnan(stats{m}.delay) && ~isinf(stats{m}.delay)
                DEL_mat(c, m) = stats{m}.delay;
            end
        end
    end

    %% ====================================================================
    %  Figure 1: Performance Metrics Grouped Bar Charts
    % =====================================================================
    fig1 = figure('Name', 'Performance Metrics Comparison', ...
        'Position', [30 30 1500 500], 'Color', 'w');

    metric_names = {'FDR (%)', 'FAR (%)', 'Accuracy (%)'};
    metric_data = {FDR_mat, FAR_mat, ACC_mat};
    % For FDR and Accuracy, higher is better; for FAR, lower is better
    metric_better = {'high', 'low', 'high'};

    for mi = 1:3
        ax = subplot(1, 3, mi);
        hold on;

        data = metric_data{mi};  % n_cond x 3
        b = bar(1:n_cond, data, 'grouped');
        for m = 1:3
            b(m).FaceColor = colors_m{m};
            b(m).EdgeColor = 'k';
            b(m).LineWidth = 0.5;
        end

        % Add value labels on bars
        for c = 1:n_cond
            for m = 1:3
                x_pos = b(m).XEndPoints(c);
                y_pos = data(c, m);
                text(x_pos, y_pos + max(data(:))*0.02, sprintf('%.1f', data(c,m)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 8, ...
                    'FontWeight', 'bold');
            end
        end

        set(gca, 'XTick', 1:n_cond, 'XTickLabel', cond_labels_short);
        ylabel(metric_names{mi}, 'FontSize', 12, 'FontWeight', 'bold');
        title(metric_names{mi}, 'FontSize', 13, 'FontWeight', 'bold');

        if strcmp(metric_better{mi}, 'high')
            ylim([0, min(110, max(data(:))*1.15)]);
        else
            ylim([0, max(data(:))*1.3 + 0.1]);
        end

        if mi == 1
            legend(methods, 'Location', 'southeast', 'FontSize', 10);
        end
        set(ax, 'FontSize', 11, 'Box', 'on');
        grid on;
        hold off;
    end

    sgtitle('Algorithm Performance Metrics Across All Conditions', ...
        'FontSize', 15, 'FontWeight', 'bold');
    saveas(fig1, 'comparison_metrics_bar.fig');
    saveas(fig1, 'comparison_metrics_bar.png');

    %% ====================================================================
    %  Figure 2: Performance Heatmap
    % =====================================================================
    fig2 = figure('Name', 'Performance Heatmap', ...
        'Position', [50 50 1100 500], 'Color', 'w');

    % Build data matrix: rows = methods, columns = metric-condition pairs
    % Columns: FDR-C1, FDR-C2, FDR-C3, FAR-C1, FAR-C2, FAR-C3, ACC-C1, ACC-C2, ACC-C3
    hm_data = [FDR_mat', FAR_mat', ACC_mat'];  % 3 x 9
    col_labels = {};
    for mi = 1:3
        mn = {'FDR','FAR','ACC'};
        for c = 1:n_cond
            col_labels{end+1} = sprintf('%s-%s', mn{mi}, cond_labels_short{c});
        end
    end

    % Custom heatmap using imagesc
    ax = axes;
    imagesc(hm_data);

    % Color map: green for good, red for bad
    n_colors = 256;
    cmap = zeros(n_colors, 3);
    for ci = 1:n_colors
        frac = (ci-1) / (n_colors-1);
        if frac < 0.5
            % Red to Yellow
            cmap(ci,:) = [0.9, 0.3+0.6*(frac/0.5), 0.2];
        else
            % Yellow to Green
            cmap(ci,:) = [0.9-0.6*((frac-0.5)/0.5), 0.9, 0.2+0.3*((frac-0.5)/0.5)];
        end
    end
    colormap(ax, cmap);
    colorbar;

    % Add text values in cells
    for r = 1:3
        for ci = 1:size(hm_data, 2)
            val = hm_data(r, ci);
            % Determine text color for readability
            if val > 50
                txt_color = 'k';
            else
                txt_color = 'w';
            end
            text(ci, r, sprintf('%.1f', val), 'HorizontalAlignment', 'center', ...
                'FontSize', 10, 'FontWeight', 'bold', 'Color', txt_color);
        end
    end

    set(gca, 'XTick', 1:length(col_labels), 'XTickLabel', col_labels, ...
        'XTickLabelRotation', 45, 'YTick', 1:3, 'YTickLabel', methods, ...
        'FontSize', 11);
    title('Performance Heatmap (% values)', 'FontSize', 14, 'FontWeight', 'bold');

    % Add separator lines between metric groups
    hold on;
    for sep = [3.5, 6.5]
        xline(sep, 'k-', 'LineWidth', 2);
    end
    hold off;

    saveas(fig2, 'comparison_heatmap.fig');
    saveas(fig2, 'comparison_heatmap.png');

    %% ====================================================================
    %  Figure 3: Radar Chart for Multi-dimensional Comparison
    % =====================================================================
    fig3 = figure('Name', 'Radar Chart Comparison', ...
        'Position', [70 70 800 750], 'Color', 'w');
    ax = axes('Position', [0.1 0.1 0.8 0.75]);

    % Metrics for radar: avg FDR, 1-avg FAR, avg Accuracy, isolation quality, detection SNR
    avg_FDR = mean(FDR_mat, 1);           % 1x3
    avg_FAR = mean(FAR_mat, 1);           % 1x3
    avg_ACC = mean(ACC_mat, 1);           % 1x3

    % Compute average detection SNR (signal-to-noise ratio in FD)
    avg_snr = zeros(1, 3);
    for m = 1:3
        snr_vals = zeros(n_cond, 1);
        for c = 1:n_cond
            FDs = {all_results(c).FD_glt, all_results(c).FD_wglt, all_results(c).FD_swglt};
            fm = all_results(c).fault_mask;
            fd_fault = FDs{m}(fm);
            fd_normal = FDs{m}(~fm);
            if std(fd_normal) > 0
                snr_vals(c) = (mean(fd_fault) - mean(fd_normal)) / std(fd_normal);
            end
        end
        avg_snr(m) = mean(snr_vals);
    end

    % Compute average isolation discrimination ratio
    avg_disc = zeros(1, 3);
    sensor_indices = [cfg.fault1.sensor_idx, cfg.fault2.sensor_idx, cfg.fault3.sensor_idx];
    for m = 1:3
        disc_vals = zeros(n_cond, 1);
        for c = 1:n_cond
            FIs = {all_results(c).FI_glt, all_results(c).FI_wglt, all_results(c).FI_swglt};
            fm = all_results(c).fault_mask;
            fi_mean = mean(FIs{m}(fm, :), 1);
            ts_idx = sensor_indices(c);
            fi_true = fi_mean(ts_idx);
            fi_others = fi_mean; fi_others(ts_idx) = 0;
            fi_max_other = max(fi_others);
            if fi_max_other > 0
                disc_vals(c) = fi_true / fi_max_other;
            else
                disc_vals(c) = 10;  % cap
            end
        end
        avg_disc(m) = min(mean(disc_vals), 10);
    end

    % Normalize all to [0, 1] range for radar chart
    radar_labels = {'FDR', '1-FAR', 'Accuracy', 'Detection SNR', 'Isolation Ratio'};
    n_axes = length(radar_labels);
    radar_data = zeros(3, n_axes);  % 3 methods x 5 metrics
    radar_data(:, 1) = avg_FDR' / 100;
    radar_data(:, 2) = (100 - avg_FAR') / 100;
    radar_data(:, 3) = avg_ACC' / 100;
    radar_data(:, 4) = min(avg_snr', max(avg_snr)) / max(avg_snr);
    radar_data(:, 5) = min(avg_disc', max(avg_disc)) / max(avg_disc);

    % Draw radar chart manually
    angles = linspace(0, 2*pi, n_axes + 1);
    angles = angles(1:end-1);

    hold on;
    % Draw grid circles
    for r = 0.2:0.2:1.0
        x_circle = r * cos(linspace(0, 2*pi, 100));
        y_circle = r * sin(linspace(0, 2*pi, 100));
        plot(x_circle, y_circle, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 0.5);
    end

    % Draw axes
    for a = 1:n_axes
        plot([0, 1.1*cos(angles(a))], [0, 1.1*sin(angles(a))], '-', ...
            'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
        text(1.2*cos(angles(a)), 1.2*sin(angles(a)), radar_labels{a}, ...
            'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end

    % Plot each method
    for m = 1:3
        r_vals = radar_data(m, :);
        x_pts = [r_vals .* cos(angles), r_vals(1) * cos(angles(1))];
        y_pts = [r_vals .* sin(angles), r_vals(1) * sin(angles(1))];
        fill(x_pts, y_pts, colors_m{m}, 'FaceAlpha', 0.15, 'EdgeColor', colors_m{m}, ...
            'LineWidth', 2.0);
        plot(r_vals .* cos(angles), r_vals .* sin(angles), 'o', ...
            'MarkerSize', 6, 'MarkerFaceColor', colors_m{m}, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end

    axis equal;
    xlim([-1.5, 1.5]); ylim([-1.5, 1.5]);
    set(gca, 'XTick', [], 'YTick', [], 'Box', 'off');
    legend(methods, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
        'FontSize', 12);
    title('Multi-dimensional Algorithm Comparison', ...
        'FontSize', 14, 'FontWeight', 'bold');
    hold off;

    saveas(fig3, 'comparison_radar.fig');
    saveas(fig3, 'comparison_radar.png');

    %% ====================================================================
    %  Figure 4: ROC Curves (FDR vs FAR at varying thresholds)
    % =====================================================================
    fig4 = figure('Name', 'ROC Curves', ...
        'Position', [90 90 1400 450], 'Color', 'w');

    for c = 1:n_cond
        ax = subplot(1, 3, c);
        hold on;

        FDs = {all_results(c).FD_glt, all_results(c).FD_wglt, all_results(c).FD_swglt};
        fm = all_results(c).fault_mask;

        for m = 1:3
            fd = FDs{m};
            fd_fault = fd(fm);
            fd_normal = fd(~fm);

            % Vary threshold from 0 to max
            n_pts = 200;
            t_range = linspace(0, max(fd)*1.2, n_pts);
            fdr_curve = zeros(n_pts, 1);
            far_curve = zeros(n_pts, 1);

            for ti = 1:n_pts
                thresh_i = t_range(ti);
                fdr_curve(ti) = sum(fd_fault > thresh_i) / max(length(fd_fault), 1);
                far_curve(ti) = sum(fd_normal > thresh_i) / max(length(fd_normal), 1);
            end

            plot(far_curve * 100, fdr_curve * 100, '-', ...
                'Color', colors_m{m}, 'LineWidth', 2.0);
        end

        % Mark operating points
        threshs = {all_results(c).T_D, all_results(c).T_D, all_results(c).T_adaptive};
        for m = 1:3
            fd = FDs{m};
            fdr_op = sum(fd(fm) > threshs{m}) / max(sum(fm), 1) * 100;
            far_op = sum(fd(~fm) > threshs{m}) / max(sum(~fm), 1) * 100;
            plot(far_op, fdr_op, 'o', 'MarkerSize', 10, ...
                'MarkerFaceColor', colors_m{m}, 'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.5);
        end

        % Ideal point
        plot(0, 100, 'p', 'MarkerSize', 14, 'MarkerFaceColor', [0 0.7 0], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1);

        xlabel('FAR (%)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('FDR (%)', 'FontSize', 11, 'FontWeight', 'bold');
        title(cond_labels{c}, 'FontSize', 12, 'FontWeight', 'bold');
        if c == 1
            legend([methods, {'Operating Point'}, {'Ideal'}], ...
                'Location', 'southeast', 'FontSize', 9);
        end
        set(ax, 'FontSize', 10, 'Box', 'on');
        grid on;
        xlim([-2, 100]); ylim([0, 105]);
        hold off;
    end

    sgtitle('ROC Curves: FDR vs FAR', 'FontSize', 15, 'FontWeight', 'bold');
    saveas(fig4, 'comparison_roc.fig');
    saveas(fig4, 'comparison_roc.png');

    %% ====================================================================
    %  Figure 5: Detection SNR and Noise Suppression Analysis
    % =====================================================================
    fig5 = figure('Name', 'Detection SNR Analysis', ...
        'Position', [110 110 1400 500], 'Color', 'w');

    % Subplot 1: Detection SNR per condition
    ax1 = subplot(1, 3, 1);
    hold on;
    snr_data = zeros(n_cond, 3);
    for c = 1:n_cond
        FDs = {all_results(c).FD_glt, all_results(c).FD_wglt, all_results(c).FD_swglt};
        fm = all_results(c).fault_mask;
        for m = 1:3
            mu_f = mean(FDs{m}(fm));
            mu_n = mean(FDs{m}(~fm));
            sig_n = std(FDs{m}(~fm));
            if sig_n > 0
                snr_data(c, m) = (mu_f - mu_n) / sig_n;
            end
        end
    end
    b = bar(1:n_cond, snr_data, 'grouped');
    for m = 1:3
        b(m).FaceColor = colors_m{m};
        b(m).EdgeColor = 'k';
    end
    % Value labels
    for c = 1:n_cond
        for m = 1:3
            text(b(m).XEndPoints(c), snr_data(c,m) + max(snr_data(:))*0.02, ...
                sprintf('%.1f', snr_data(c,m)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
    set(gca, 'XTick', 1:n_cond, 'XTickLabel', cond_labels_short);
    ylabel('SNR (dB-like)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Detection SNR', 'FontSize', 13, 'FontWeight', 'bold');
    legend(methods, 'Location', 'best', 'FontSize', 10);
    set(ax1, 'FontSize', 11, 'Box', 'on');
    grid on; hold off;

    % Subplot 2: Normal-period FD variance comparison
    ax2 = subplot(1, 3, 2);
    hold on;
    var_data = zeros(n_cond, 3);
    for c = 1:n_cond
        FDs = {all_results(c).FD_glt, all_results(c).FD_wglt, all_results(c).FD_swglt};
        fm = all_results(c).fault_mask;
        for m = 1:3
            var_data(c, m) = var(FDs{m}(~fm));
        end
    end
    b2 = bar(1:n_cond, var_data, 'grouped');
    for m = 1:3
        b2(m).FaceColor = colors_m{m};
        b2(m).EdgeColor = 'k';
    end
    for c = 1:n_cond
        for m = 1:3
            text(b2(m).XEndPoints(c), var_data(c,m) + max(var_data(:))*0.02, ...
                sprintf('%.1f', var_data(c,m)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
    set(gca, 'XTick', 1:n_cond, 'XTickLabel', cond_labels_short);
    ylabel('Var(FD) normal', 'FontSize', 11, 'FontWeight', 'bold');
    title('Normal-Period FD Variance', 'FontSize', 13, 'FontWeight', 'bold');
    legend(methods, 'Location', 'best', 'FontSize', 10);
    set(ax2, 'FontSize', 11, 'Box', 'on');
    grid on; hold off;

    % Subplot 3: Fault-period mean FD / threshold ratio
    ax3 = subplot(1, 3, 3);
    hold on;
    ratio_data = zeros(n_cond, 3);
    for c = 1:n_cond
        FDs = {all_results(c).FD_glt, all_results(c).FD_wglt, all_results(c).FD_swglt};
        threshs = {all_results(c).T_D, all_results(c).T_D, all_results(c).T_adaptive};
        fm = all_results(c).fault_mask;
        for m = 1:3
            ratio_data(c, m) = mean(FDs{m}(fm)) / threshs{m};
        end
    end
    b3 = bar(1:n_cond, ratio_data, 'grouped');
    for m = 1:3
        b3(m).FaceColor = colors_m{m};
        b3(m).EdgeColor = 'k';
    end
    for c = 1:n_cond
        for m = 1:3
            text(b3(m).XEndPoints(c), ratio_data(c,m) + max(ratio_data(:))*0.02, ...
                sprintf('%.1f', ratio_data(c,m)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
    yline(1, 'k--', 'LineWidth', 1.5);
    set(gca, 'XTick', 1:n_cond, 'XTickLabel', cond_labels_short);
    ylabel('Mean FD_{fault} / Threshold', 'FontSize', 11, 'FontWeight', 'bold');
    title('Threshold Exceedance Ratio', 'FontSize', 13, 'FontWeight', 'bold');
    legend(methods, 'Location', 'best', 'FontSize', 10);
    set(ax3, 'FontSize', 11, 'Box', 'on');
    grid on; hold off;

    sgtitle('Detection Quality Analysis', 'FontSize', 15, 'FontWeight', 'bold');
    saveas(fig5, 'comparison_snr.fig');
    saveas(fig5, 'comparison_snr.png');

    %% ====================================================================
    %  Figure 6: Comprehensive Summary Dashboard
    % =====================================================================
    fig6 = figure('Name', 'Summary Dashboard', ...
        'Position', [20 20 1600 1000], 'Color', 'w');

    % --- Panel 1: Performance table (top-left) ---
    ax1 = subplot(2, 3, 1);
    axis off;

    % Build table data
    tbl_headers = {'', 'FDR(%)', 'FAR(%)', 'ACC(%)'};
    tbl_rows = {};
    row_colors = {};
    for c = 1:n_cond
        for m = 1:3
            tbl_rows{end+1} = {sprintf('%s-%s', cond_labels_short{c}, methods{m}), ...
                sprintf('%.1f', FDR_mat(c,m)), ...
                sprintf('%.2f', FAR_mat(c,m)), ...
                sprintf('%.1f', ACC_mat(c,m))};
            row_colors{end+1} = colors_m{m};
        end
    end

    % Draw table
    n_rows = length(tbl_rows);
    n_cols = 4;
    row_h = 1 / (n_rows + 2);
    col_w = 1 / n_cols;

    % Headers
    for j = 1:n_cols
        text(col_w*(j-0.5), 1 - row_h*0.5, tbl_headers{j}, ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', 9, 'FontWeight', 'bold');
    end
    % Draw header line
    line([0 1], [1-row_h, 1-row_h], 'Color', 'k', 'LineWidth', 1, ...
        'Clipping', 'off');

    for ri = 1:n_rows
        y_pos = 1 - row_h * (ri + 0.5);
        for j = 1:n_cols
            tc = 'k';
            fw = 'normal';
            if j == 1
                tc = row_colors{ri};
                fw = 'bold';
            end
            text(col_w*(j-0.5), y_pos, tbl_rows{ri}{j}, ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontSize', 8, 'Color', tc, 'FontWeight', fw);
        end
    end
    title('Performance Summary', 'FontSize', 12, 'FontWeight', 'bold');

    % --- Panel 2: Average performance bar (top-center) ---
    ax2 = subplot(2, 3, 2);
    hold on;
    avg_data = [mean(FDR_mat, 1)', mean(100-FAR_mat, 1)', mean(ACC_mat, 1)'];
    b = bar(1:3, avg_data, 'grouped');
    bar_colors = {[0.3 0.7 0.3], [0.3 0.5 0.8], [0.8 0.6 0.2]};
    for j = 1:3
        b(j).FaceColor = bar_colors{j};
        b(j).EdgeColor = 'k';
    end
    for mi = 1:3
        for j = 1:3
            text(b(j).XEndPoints(mi), avg_data(mi,j) + 1, ...
                sprintf('%.1f', avg_data(mi,j)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
    set(gca, 'XTick', 1:3, 'XTickLabel', methods);
    ylabel('%', 'FontSize', 10);
    title('Average Performance', 'FontSize', 12, 'FontWeight', 'bold');
    legend({'Avg FDR', 'Avg (1-FAR)', 'Avg Accuracy'}, ...
        'Location', 'southeast', 'FontSize', 8);
    ylim([0, 110]);
    set(ax2, 'FontSize', 10, 'Box', 'on');
    grid on; hold off;

    % --- Panel 3: Improvement over GLT (top-right) ---
    ax3 = subplot(2, 3, 3);
    hold on;
    % Compute improvement of WGLT and SWGLT over GLT
    improv_labels = {'FDR gain', 'FAR reduction', 'ACC gain'};
    improv_wglt = [mean(FDR_mat(:,2) - FDR_mat(:,1)), ...
                   mean(FAR_mat(:,1) - FAR_mat(:,2)), ...
                   mean(ACC_mat(:,2) - ACC_mat(:,1))];
    improv_swglt = [mean(FDR_mat(:,3) - FDR_mat(:,1)), ...
                    mean(FAR_mat(:,1) - FAR_mat(:,3)), ...
                    mean(ACC_mat(:,3) - ACC_mat(:,1))];

    x = 1:3;
    bar_w = 0.35;
    b1 = bar(x - bar_w/2, improv_wglt, bar_w, 'FaceColor', color_wglt, 'EdgeColor', 'k');
    b2 = bar(x + bar_w/2, improv_swglt, bar_w, 'FaceColor', color_swglt, 'EdgeColor', 'k');

    % Value labels
    for j = 1:3
        text(x(j) - bar_w/2, improv_wglt(j) + sign(improv_wglt(j))*0.5, ...
            sprintf('%+.1f', improv_wglt(j)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
        text(x(j) + bar_w/2, improv_swglt(j) + sign(improv_swglt(j))*0.5, ...
            sprintf('%+.1f', improv_swglt(j)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
    end

    yline(0, 'k-', 'LineWidth', 1);
    set(gca, 'XTick', x, 'XTickLabel', improv_labels);
    ylabel('Percentage Points', 'FontSize', 10);
    title('Improvement over GLT', 'FontSize', 12, 'FontWeight', 'bold');
    legend({'WGLT', 'SWGLT'}, 'Location', 'best', 'FontSize', 10);
    set(ax3, 'FontSize', 10, 'Box', 'on');
    grid on; hold off;

    % --- Panel 4: Mini radar (bottom-left) ---
    ax4 = subplot(2, 3, 4);
    hold on;
    % Reuse radar data computed earlier
    radar_data_local = zeros(3, n_axes);
    radar_data_local(:, 1) = avg_FDR' / 100;
    radar_data_local(:, 2) = (100 - avg_FAR') / 100;
    radar_data_local(:, 3) = avg_ACC' / 100;
    radar_data_local(:, 4) = min(avg_snr', max(avg_snr)) / max(max(avg_snr), 1e-10);
    radar_data_local(:, 5) = min(avg_disc', max(avg_disc)) / max(max(avg_disc), 1e-10);

    for r = 0.25:0.25:1.0
        x_c = r * cos(linspace(0, 2*pi, 100));
        y_c = r * sin(linspace(0, 2*pi, 100));
        plot(x_c, y_c, '-', 'Color', [0.88 0.88 0.88], 'LineWidth', 0.4);
    end
    for a = 1:n_axes
        plot([0, cos(angles(a))], [0, sin(angles(a))], '-', ...
            'Color', [0.8 0.8 0.8], 'LineWidth', 0.5);
        text(1.15*cos(angles(a)), 1.15*sin(angles(a)), radar_labels{a}, ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
    end
    for m = 1:3
        r_vals = radar_data_local(m, :);
        x_pts = [r_vals .* cos(angles), r_vals(1)*cos(angles(1))];
        y_pts = [r_vals .* sin(angles), r_vals(1)*sin(angles(1))];
        fill(x_pts, y_pts, colors_m{m}, 'FaceAlpha', 0.12, ...
            'EdgeColor', colors_m{m}, 'LineWidth', 1.5);
    end
    axis equal; xlim([-1.4, 1.4]); ylim([-1.4, 1.4]);
    set(gca, 'XTick', [], 'YTick', []);
    title('Capability Radar', 'FontSize', 12, 'FontWeight', 'bold');
    hold off;

    % --- Panel 5: FDR vs FAR scatter (bottom-center) ---
    ax5 = subplot(2, 3, 5);
    hold on;
    marker_shapes = {'o', 's', 'd'};  % circle, square, diamond for conditions
    for c = 1:n_cond
        for m = 1:3
            plot(FAR_mat(c,m), FDR_mat(c,m), marker_shapes{c}, ...
                'MarkerSize', 12, 'MarkerFaceColor', colors_m{m}, ...
                'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
        end
    end
    % Ideal region
    patch([0 2 2 0], [95 95 105 105], [0.9 1 0.9], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.4);

    xlabel('FAR (%)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('FDR (%)', 'FontSize', 11, 'FontWeight', 'bold');
    title('FDR vs FAR Operating Points', 'FontSize', 12, 'FontWeight', 'bold');

    % Create custom legend entries
    h_leg = [];
    for m = 1:3
        h_leg(end+1) = plot(NaN, NaN, 'o', 'MarkerSize', 8, ...
            'MarkerFaceColor', colors_m{m}, 'MarkerEdgeColor', 'k');
    end
    for c = 1:n_cond
        h_leg(end+1) = plot(NaN, NaN, marker_shapes{c}, 'MarkerSize', 8, ...
            'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k');
    end
    legend(h_leg, [methods, cond_labels_short], 'Location', 'best', 'FontSize', 8);

    set(ax5, 'FontSize', 10, 'Box', 'on');
    grid on; hold off;

    % --- Panel 6: Detection delay for soft fault (bottom-right) ---
    ax6 = subplot(2, 3, 6);
    hold on;
    % Only condition 3 has meaningful delay
    if any(~isnan(DEL_mat(3, :)))
        delays = DEL_mat(3, :);
        delays(isnan(delays)) = 0;
        delays(isinf(delays)) = max(delays(~isinf(delays))) * 1.5;

        b_del = bar(1:3, delays);
        for m = 1:3
            b_del.FaceColor = 'flat';
            b_del.CData(m, :) = colors_m{m};
        end
        b_del.EdgeColor = 'k';

        for m = 1:3
            if ~isnan(DEL_mat(3, m)) && ~isinf(DEL_mat(3, m))
                text(m, delays(m) + max(delays)*0.05, sprintf('%.1fs', delays(m)), ...
                    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
            else
                text(m, delays(m) + max(delays)*0.05, 'N/A', ...
                    'HorizontalAlignment', 'center', 'FontSize', 10, ...
                    'FontWeight', 'bold', 'Color', [0.8 0 0]);
            end
        end
        set(gca, 'XTick', 1:3, 'XTickLabel', methods);
        ylabel('Delay (s)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Soft Fault Detection Delay (Cond 3)', 'FontSize', 12, 'FontWeight', 'bold');
    else
        text(0.5, 0.5, 'No delay data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 12);
        title('Detection Delay', 'FontSize', 12, 'FontWeight', 'bold');
    end
    set(ax6, 'FontSize', 10, 'Box', 'on');
    grid on; hold off;

    sgtitle('Comprehensive Algorithm Comparison Dashboard', ...
        'FontSize', 16, 'FontWeight', 'bold');
    saveas(fig6, 'comparison_dashboard.fig');
    saveas(fig6, 'comparison_dashboard.png');

    %% Print summary to console
    fprintf('\n====================================================================\n');
    fprintf('  COMPREHENSIVE COMPARISON SUMMARY\n');
    fprintf('====================================================================\n');
    fprintf('\n  Average Performance (across all conditions):\n');
    fprintf('  %-8s  FDR=%6.1f%%  FAR=%6.2f%%  Accuracy=%6.1f%%\n', ...
        'GLT',   mean(FDR_mat(:,1)), mean(FAR_mat(:,1)), mean(ACC_mat(:,1)));
    fprintf('  %-8s  FDR=%6.1f%%  FAR=%6.2f%%  Accuracy=%6.1f%%\n', ...
        'WGLT',  mean(FDR_mat(:,2)), mean(FAR_mat(:,2)), mean(ACC_mat(:,2)));
    fprintf('  %-8s  FDR=%6.1f%%  FAR=%6.2f%%  Accuracy=%6.1f%%\n', ...
        'SWGLT', mean(FDR_mat(:,3)), mean(FAR_mat(:,3)), mean(ACC_mat(:,3)));

    fprintf('\n  SWGLT Improvement over GLT:\n');
    fprintf('    FDR:      %+.1f pp\n', mean(FDR_mat(:,3) - FDR_mat(:,1)));
    fprintf('    FAR:      %.2f%% -> %.2f%% (%.1fx reduction)\n', ...
        mean(FAR_mat(:,1)), mean(FAR_mat(:,3)), ...
        mean(FAR_mat(:,1)) / max(mean(FAR_mat(:,3)), 0.001));
    fprintf('    Accuracy: %+.1f pp\n', mean(ACC_mat(:,3) - ACC_mat(:,1)));
    fprintf('====================================================================\n');

end
