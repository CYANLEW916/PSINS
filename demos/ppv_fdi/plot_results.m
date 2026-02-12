function plot_results(t, FD_glt, FD_wglt, FD_swglt, ...
    T_D, T_adaptive, ...
    FI_glt, FI_wglt, FI_swglt, ...
    fault_mask, fault_intervals, ...
    cond_num, cond_title, cfg)
% PLOT_RESULTS  Plot fault detection and isolation results.
%   PLOT_RESULTS(...) generates detection function time series and isolation
%   bar charts for each fault condition comparing GLT, WGLT, and SWGLT.
%
% See also  main_simulation, evaluate_performance.

    sensor_labels = {'INS1-X','INS1-Y','INS1-Z', ...
                     'INS2-X','INS2-Y','INS2-Z', ...
                     'ISIS-X','ISIS-Y','ISIS-Z'};

    colors = struct('glt', [0.2 0.4 0.8], ...
                    'wglt', [0.8 0.5 0.1], ...
                    'swglt', [0.1 0.7 0.3]);

    %% Figure 1: Detection function time series (3 subplots)
    fig1 = figure('Name', sprintf('Condition %d - Detection Functions', cond_num), ...
        'Position', [50 50 1200 800]);

    methods = {'GLT', 'WGLT', 'SWGLT'};
    FDs = {FD_glt, FD_wglt, FD_swglt};
    thresholds = {T_D, T_D, T_adaptive};
    method_colors = {colors.glt, colors.wglt, colors.swglt};

    for m = 1:3
        subplot(3, 1, m);
        hold on;

        % Shade fault intervals
        yl = [0, max(FDs{m}) * 1.2];
        if yl(2) <= yl(1), yl(2) = yl(1) + 1; end
        for fi = 1:size(fault_intervals, 1)
            patch([fault_intervals(fi,1) fault_intervals(fi,2) ...
                   fault_intervals(fi,2) fault_intervals(fi,1)], ...
                  [yl(1) yl(1) yl(2) yl(2)], ...
                  [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end

        % Plot detection function
        plot(t, FDs{m}, 'Color', method_colors{m}, 'LineWidth', 0.5);

        % Plot threshold
        if isscalar(thresholds{m})
            yline(thresholds{m}, 'r--', 'LineWidth', 1.5);
        else
            plot(t, ones(size(t)) * thresholds{m}, 'r--', 'LineWidth', 1.5);
        end

        ylabel('FD');
        title(sprintf('%s - %s', methods{m}, cond_title));
        xlim([t(1) t(end)]);
        ylim(yl);
        if m == 3
            xlabel('Time (s)');
        end
        legend('Fault Period', 'Detection Function', 'Threshold', ...
            'Location', 'northeast');
        grid on;
        hold off;
    end

    saveas(fig1, fullfile(fileparts(mfilename('fullpath')), ...
        sprintf('cond%d_detection.fig', cond_num)));
    saveas(fig1, fullfile(fileparts(mfilename('fullpath')), ...
        sprintf('cond%d_detection.png', cond_num)));

    %% Figure 2: Fault isolation results
    fig2 = figure('Name', sprintf('Condition %d - Fault Isolation', cond_num), ...
        'Position', [100 100 1200 600]);

    % For isolation, take the mean FI during fault periods
    fault_idx = find(fault_mask);
    if isempty(fault_idx)
        return;
    end

    FIs = {FI_glt, FI_wglt, FI_swglt};

    for m = 1:3
        subplot(1, 3, m);

        % Average isolation metric during fault periods
        FI_fault_mean = mean(FIs{m}(fault_idx, :), 1);

        bar(1:9, FI_fault_mean, 'FaceColor', method_colors{m}, 'EdgeColor', 'k');
        set(gca, 'XTick', 1:9, 'XTickLabel', sensor_labels, 'XTickLabelRotation', 45);
        ylabel('FI (mean during fault)');
        title(sprintf('%s Isolation', methods{m}));
        grid on;

        % Highlight the max (isolated sensor)
        [~, max_idx] = max(FI_fault_mean);
        hold on;
        bar(max_idx, FI_fault_mean(max_idx), 'FaceColor', [1 0.3 0.3], 'EdgeColor', 'k');
        hold off;
    end

    sgtitle(cond_title);
    saveas(fig2, fullfile(fileparts(mfilename('fullpath')), ...
        sprintf('cond%d_isolation.fig', cond_num)));
    saveas(fig2, fullfile(fileparts(mfilename('fullpath')), ...
        sprintf('cond%d_isolation.png', cond_num)));

end
