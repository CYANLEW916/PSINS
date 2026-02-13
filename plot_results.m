function plot_results(t, FD_glt, FD_wglt, FD_swglt, ...
    T_D, T_adaptive, ...
    FI_glt, FI_wglt, FI_swglt, ...
    fault_mask, fault_intervals, ...
    cond_num, cond_title, cfg)
% PLOT_RESULTS  Comprehensive per-condition FDI visualization.
%   Generates four figures per fault condition:
%     Figure 1: Enhanced detection function time series (GLT/WGLT/SWGLT)
%     Figure 2: Zoomed detection around each fault interval
%     Figure 3: False alarm histogram analysis
%     Figure 4: Grouped isolation comparison bar chart
%
% See also  main_simulation, evaluate_performance, plot_comprehensive_comparison.

    sensor_labels = {'INS1-X','INS1-Y','INS1-Z', ...
                     'INS2-X','INS2-Y','INS2-Z', ...
                     'ISIS-X','ISIS-Y','ISIS-Z'};

    % Color scheme: professional and distinguishable
    color_glt  = [0.55 0.55 0.55];   % gray for traditional GLT
    color_wglt = [0.20 0.47 0.87];   % blue for WGLT
    color_swglt= [0.85 0.22 0.16];   % red for SWGLT (highlight)

    colors_m = {color_glt, color_wglt, color_swglt};
    methods = {'GLT', 'WGLT', 'SWGLT'};
    FDs = {FD_glt, FD_wglt, FD_swglt};
    thresholds = {T_D, T_D, T_adaptive};
    thresh_labels = {'T_D (fixed)', 'T_D (fixed)', 'T_{adaptive}'};

    % Determine true faulty sensor from config
    switch cond_num
        case 1, true_sensor = cfg.fault1.sensor_idx;
        case 2, true_sensor = cfg.fault2.sensor_idx;
        case 3, true_sensor = cfg.fault3.sensor_idx;
    end

    %% ====================================================================
    %  Figure 1: Enhanced Detection Function Time Series
    % =====================================================================
    fig1 = figure('Name', sprintf('Condition %d - Detection', cond_num), ...
        'Position', [30 30 1400 900], 'Color', 'w');

    for m = 1:3
        ax = subplot(3, 1, m);
        hold on;

        % Y-axis range
        fd_max = max(FDs{m});
        yl = [0, max([fd_max, thresholds{m}]) * 1.4];
        if yl(2) <= yl(1), yl(2) = yl(1) + 1; end

        % Shade fault intervals (light red)
        for fi = 1:size(fault_intervals, 1)
            patch([fault_intervals(fi,1) fault_intervals(fi,2) ...
                   fault_intervals(fi,2) fault_intervals(fi,1)], ...
                  [yl(1) yl(1) yl(2) yl(2)], ...
                  [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end

        % Detection function
        plot(t, FDs{m}, 'Color', colors_m{m}, 'LineWidth', 1.0);

        % Threshold line
        yline(thresholds{m}, 'k--', 'LineWidth', 1.8);

        % Add threshold value text
        text(t(end)*0.98, thresholds{m}*1.08, ...
            sprintf('%s = %.2f', thresh_labels{m}, thresholds{m}), ...
            'HorizontalAlignment', 'right', 'FontSize', 9, ...
            'FontWeight', 'bold', 'Color', [0.3 0.3 0.3]);

        ylabel('FD', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s', methods{m}), 'FontSize', 13, 'FontWeight', 'bold');
        xlim([t(1) t(end)]);
        ylim(yl);
        if m == 3
            xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
        end
        legend('Fault Period', 'FD(t)', 'Threshold', ...
            'Location', 'northeast', 'FontSize', 9);
        set(ax, 'FontSize', 10, 'LineWidth', 0.8, 'Box', 'on');
        grid on; grid minor;
        hold off;
    end

    sgtitle(cond_title, 'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig1, sprintf('cond%d_detection.fig', cond_num));
    saveas(fig1, sprintf('cond%d_detection.png', cond_num));

    %% ====================================================================
    %  Figure 2: Zoomed Detection Around Fault Intervals
    % =====================================================================
    n_intervals = size(fault_intervals, 1);
    fig2 = figure('Name', sprintf('Condition %d - Zoomed Detection', cond_num), ...
        'Position', [50 50 1500 300*n_intervals], 'Color', 'w');

    margin_sec = 10;  % seconds of context before/after fault

    for fi = 1:n_intervals
        t_start = fault_intervals(fi, 1);
        t_end   = fault_intervals(fi, 2);
        zoom_mask = (t >= t_start - margin_sec) & (t <= t_end + margin_sec);
        t_zoom = t(zoom_mask);

        for m = 1:3
            ax = subplot(n_intervals, 3, (fi-1)*3 + m);
            hold on;

            fd_zoom = FDs{m}(zoom_mask);
            yl_z = [0, max([max(fd_zoom), thresholds{m}]) * 1.5];
            if yl_z(2) <= yl_z(1), yl_z(2) = yl_z(1) + 1; end

            % Fault interval shading
            patch([t_start t_end t_end t_start], ...
                  [yl_z(1) yl_z(1) yl_z(2) yl_z(2)], ...
                  [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

            % Detection function
            plot(t_zoom, fd_zoom, 'Color', colors_m{m}, 'LineWidth', 1.2);

            % Threshold
            yline(thresholds{m}, 'k--', 'LineWidth', 1.5);

            % Mark first detection point within fault
            fault_zoom_mask = zoom_mask & fault_mask;
            detected_in_fault = fault_zoom_mask & (FDs{m} > thresholds{m});
            first_det = find(detected_in_fault, 1, 'first');
            if ~isempty(first_det)
                plot(t(first_det), FDs{m}(first_det), 'v', ...
                    'MarkerSize', 8, 'MarkerFaceColor', [0 0.7 0], ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
            end

            ylim(yl_z);
            set(ax, 'FontSize', 9, 'Box', 'on');
            grid on;

            if fi == 1
                title(methods{m}, 'FontSize', 12, 'FontWeight', 'bold');
            end
            if m == 1
                ylabel(sprintf('Interval %d\nFD', fi), 'FontSize', 10, 'FontWeight', 'bold');
            end
            if fi == n_intervals
                xlabel('Time (s)', 'FontSize', 10);
            end
            hold off;
        end
    end

    sgtitle(sprintf('%s - Zoomed View', cond_title), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig2, sprintf('cond%d_detection_zoom.fig', cond_num));
    saveas(fig2, sprintf('cond%d_detection_zoom.png', cond_num));

    %% ====================================================================
    %  Figure 3: False Alarm Histogram Analysis
    % =====================================================================
    fig3 = figure('Name', sprintf('Condition %d - False Alarm Analysis', cond_num), ...
        'Position', [70 70 1400 450], 'Color', 'w');

    normal_mask = ~fault_mask;

    for m = 1:3
        ax = subplot(1, 3, m);
        hold on;

        fd_normal = FDs{m}(normal_mask);

        % Histogram of FD values during normal period
        n_bins = 80;
        [counts, edges] = histcounts(fd_normal, n_bins);
        bin_centers = (edges(1:end-1) + edges(2:end)) / 2;
        bar(bin_centers, counts, 1.0, 'FaceColor', colors_m{m}, ...
            'FaceAlpha', 0.7, 'EdgeColor', 'none');

        % Threshold line
        xline(thresholds{m}, 'k--', 'LineWidth', 2);

        % Mark false alarm region
        fa_count = sum(fd_normal > thresholds{m});
        fa_rate = fa_count / length(fd_normal) * 100;

        % Shade false alarm region
        fa_bins = bin_centers > thresholds{m};
        if any(fa_bins)
            bar(bin_centers(fa_bins), counts(fa_bins), 1.0, ...
                'FaceColor', [1 0.3 0.3], 'FaceAlpha', 0.9, 'EdgeColor', 'none');
        end

        % Annotation
        text(0.95, 0.90, sprintf('FAR = %.2f%%', fa_rate), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
        text(0.95, 0.80, sprintf('N_{FA} = %d', fa_count), ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'FontSize', 10, 'Color', [0.4 0.4 0.4]);

        title(methods{m}, 'FontSize', 13, 'FontWeight', 'bold');
        xlabel('FD value', 'FontSize', 10);
        if m == 1
            ylabel('Count', 'FontSize', 11, 'FontWeight', 'bold');
        end
        set(ax, 'FontSize', 10, 'Box', 'on');
        grid on;
        hold off;
    end

    sgtitle(sprintf('%s - Normal Period FD Distribution', cond_title), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig3, sprintf('cond%d_false_alarm.fig', cond_num));
    saveas(fig3, sprintf('cond%d_false_alarm.png', cond_num));

    %% ====================================================================
    %  Figure 4: Grouped Isolation Comparison
    % =====================================================================
    fault_idx = find(fault_mask);
    if isempty(fault_idx)
        return;
    end

    fig4 = figure('Name', sprintf('Condition %d - Isolation', cond_num), ...
        'Position', [90 90 1200 550], 'Color', 'w');

    FIs = {FI_glt, FI_wglt, FI_swglt};

    % Compute mean FI during fault periods for each method
    FI_means = zeros(3, 9);
    for m = 1:3
        FI_means(m, :) = mean(FIs{m}(fault_idx, :), 1);
    end

    % Normalize each method's FI by its max for better comparison
    FI_norm = zeros(3, 9);
    for m = 1:3
        mx = max(FI_means(m, :));
        if mx > 0
            FI_norm(m, :) = FI_means(m, :) / mx;
        end
    end

    % --- Subplot 1: Grouped bar chart (absolute values) ---
    ax1 = subplot(1, 2, 1);
    hold on;

    bar_data = FI_means';  % 9 x 3
    b = bar(1:9, bar_data, 'grouped');
    for m = 1:3
        b(m).FaceColor = colors_m{m};
        b(m).EdgeColor = 'k';
        b(m).LineWidth = 0.5;
    end

    % Highlight true faulty sensor
    xline(true_sensor, 'r-', 'LineWidth', 2, 'Alpha', 0.5);
    text(true_sensor, max(FI_means(:))*1.05, '\downarrow True Fault', ...
        'HorizontalAlignment', 'center', 'FontSize', 10, ...
        'FontWeight', 'bold', 'Color', [0.8 0 0]);

    set(gca, 'XTick', 1:9, 'XTickLabel', sensor_labels, 'XTickLabelRotation', 45);
    ylabel('FI (mean)', 'FontSize', 11, 'FontWeight', 'bold');
    title('Absolute Isolation Values', 'FontSize', 12, 'FontWeight', 'bold');
    legend(methods, 'Location', 'best', 'FontSize', 10);
    set(ax1, 'FontSize', 10, 'Box', 'on');
    grid on;
    hold off;

    % --- Subplot 2: Normalized isolation + discrimination ratio ---
    ax2 = subplot(1, 2, 2);
    hold on;

    bar_data_norm = FI_norm';  % 9 x 3
    b2 = bar(1:9, bar_data_norm, 'grouped');
    for m = 1:3
        b2(m).FaceColor = colors_m{m};
        b2(m).EdgeColor = 'k';
        b2(m).LineWidth = 0.5;
    end

    % Highlight true faulty sensor
    xline(true_sensor, 'r-', 'LineWidth', 2, 'Alpha', 0.5);

    % Compute discrimination ratio: FI(true) / max(FI(others))
    disc_text = '';
    for m = 1:3
        fi_true = FI_means(m, true_sensor);
        fi_others = FI_means(m, :);
        fi_others(true_sensor) = 0;
        fi_max_other = max(fi_others);
        if fi_max_other > 0
            disc_ratio = fi_true / fi_max_other;
        else
            disc_ratio = Inf;
        end
        disc_text = [disc_text, sprintf('%s: %.1f  ', methods{m}, disc_ratio)];
    end
    text(0.5, 0.95, ['Discrimination Ratio: ' disc_text], ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.5]);

    set(gca, 'XTick', 1:9, 'XTickLabel', sensor_labels, 'XTickLabelRotation', 45);
    ylabel('Normalized FI', 'FontSize', 11, 'FontWeight', 'bold');
    title('Normalized Isolation (each method / max)', 'FontSize', 12, 'FontWeight', 'bold');
    legend(methods, 'Location', 'best', 'FontSize', 10);
    set(ax2, 'FontSize', 10, 'Box', 'on');
    grid on;
    hold off;

    sgtitle(sprintf('%s - Fault Isolation', cond_title), ...
        'FontSize', 14, 'FontWeight', 'bold');
    saveas(fig4, sprintf('cond%d_isolation.fig', cond_num));
    saveas(fig4, sprintf('cond%d_isolation.png', cond_num));

end
