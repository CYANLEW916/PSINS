function plot_results(t, FD_glt, FD_wglt, FD_swglt, ...
    T_D, T_adaptive, ...
    FI_glt, FI_wglt, FI_swglt, ...
    fault_mask, fault_intervals, ...
    cond_num, cond_title, cfg)
% PLOT_RESULTS  Plot fault detection and isolation results.
%   Generates two figures per fault condition:
%     Figure 1: Detection function time series (GLT/WGLT/SWGLT) with threshold
%               and fault-period shading
%     Figure 2: Fault isolation bar charts (mean FI during fault periods)
%
% See also  main_simulation, evaluate_performance.

    sensor_labels = {'INS1-X','INS1-Y','INS1-Z', ...
                     'INS2-X','INS2-Y','INS2-Z', ...
                     'ISIS-X','ISIS-Y','ISIS-Z'};

    colors_m = {[0.2 0.4 0.8], [0.8 0.5 0.1], [0.1 0.7 0.3]};
    methods = {'GLT', 'WGLT', 'SWGLT'};
    FDs = {FD_glt, FD_wglt, FD_swglt};
    thresholds = {T_D, T_D, T_adaptive};

    %% Figure 1: Detection function time series
    figure('Name', sprintf('Condition %d - Detection', cond_num), ...
        'Position', [50 50 1200 800]);

    for m = 1:3
        subplot(3, 1, m);
        hold on;

        % Y-axis range
        yl = [0, max([max(FDs{m}), thresholds{m}]) * 1.3];
        if yl(2) <= yl(1), yl(2) = yl(1) + 1; end

        % Shade fault intervals (light red)
        for fi = 1:size(fault_intervals, 1)
            patch([fault_intervals(fi,1) fault_intervals(fi,2) ...
                   fault_intervals(fi,2) fault_intervals(fi,1)], ...
                  [yl(1) yl(1) yl(2) yl(2)], ...
                  [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end

        % Detection function
        plot(t, FDs{m}, 'Color', colors_m{m}, 'LineWidth', 0.5);

        % Threshold (red dashed)
        yline(thresholds{m}, 'r--', 'LineWidth', 1.5);

        ylabel('FD');
        title(sprintf('%s - %s', methods{m}, cond_title));
        xlim([t(1) t(end)]);
        ylim(yl);
        if m == 3
            xlabel('Time (s)');
        end
        legend('Fault Period', 'Detection Function', 'Threshold', 'Location', 'best');
        grid on;
        hold off;
    end

    saveas(gcf, sprintf('cond%d_detection.fig', cond_num));
    saveas(gcf, sprintf('cond%d_detection.png', cond_num));

    %% Figure 2: Fault isolation bar charts
    fault_idx = find(fault_mask);
    if isempty(fault_idx)
        return;
    end

    figure('Name', sprintf('Condition %d - Isolation', cond_num), ...
        'Position', [100 100 1200 500]);

    FIs = {FI_glt, FI_wglt, FI_swglt};

    for m = 1:3
        subplot(1, 3, m);

        FI_mean = mean(FIs{m}(fault_idx, :), 1);

        b = bar(1:9, FI_mean, 'FaceColor', colors_m{m}, 'EdgeColor', 'k');
        set(gca, 'XTick', 1:9, 'XTickLabel', sensor_labels, 'XTickLabelRotation', 45);
        ylabel('FI (mean during fault)');
        title(sprintf('%s Isolation', methods{m}));
        grid on;

        % Highlight isolated sensor (max)
        [~, max_idx] = max(FI_mean);
        hold on;
        bar(max_idx, FI_mean(max_idx), 'FaceColor', [1 0.3 0.3], 'EdgeColor', 'k');
        hold off;
    end

    sgtitle(cond_title);
    saveas(gcf, sprintf('cond%d_isolation.fig', cond_num));
    saveas(gcf, sprintf('cond%d_isolation.png', cond_num));

end
