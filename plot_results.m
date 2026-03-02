function plot_results(t, FD_glt, FD_wglt, FD_swglt, ...
    T_D, T_adaptive, ...
    FI_glt, FI_wglt, FI_swglt, ...
    fault_mask, fault_intervals, ...
    cond_num, cond_title, cfg, results_swglt, true_fault_sensor)
% PLOT_RESULTS  Plot fault detection, isolation, and SWGLT diagnostics.
%   Generates:
%     (1) Detection curves for GLT/WGLT/SWGLT
%     (2) Isolation bar charts using FI means
%     (3) SWGLT isolation robustness diagnostics (margin + posterior)

    sensor_labels = {'INS1-X','INS1-Y','INS1-Z', ...
                     'INS2-X','INS2-Y','INS2-Z', ...
                     'ISIS-X','ISIS-Y','ISIS-Z'};

    colors_m = {[0.2 0.4 0.8], [0.8 0.5 0.1], [0.1 0.7 0.3]};
    methods = {'GLT', 'WGLT', 'SWGLT'};
    FDs = {FD_glt, FD_wglt, FD_swglt};
    thresholds = {T_D, T_D, T_adaptive};

    figure('Name', sprintf('Condition %d - Detection', cond_num), ...
        'Position', [50 50 1200 800]);

    for m = 1:3
        subplot(3, 1, m);
        hold on;
        hFault = [];

        yl = [0, max([max(FDs{m}), thresholds{m}]) * 1.3];
        if yl(2) <= yl(1)
            yl(2) = yl(1) + 1;
        end

        for fi = 1:size(fault_intervals, 1)
            hPatch = patch([fault_intervals(fi,1) fault_intervals(fi,2) ...
                            fault_intervals(fi,2) fault_intervals(fi,1)], ...
                           [yl(1) yl(1) yl(2) yl(2)], ...
                           [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            if fi == 1
                hFault = hPatch;
            else
                set(hPatch, 'HandleVisibility', 'off');
            end
        end

        hFD = plot(t, FDs{m}, 'Color', colors_m{m}, 'LineWidth', 0.5);
        hTh = yline(thresholds{m}, 'r--', 'LineWidth', 1.5);

        ylabel('FD');
        if strcmp(methods{m}, 'GLT')
            title(sprintf('%s - %s (uniform sigma assumption)', methods{m}, cond_title));
        else
            title(sprintf('%s - %s', methods{m}, cond_title));
        end

        xlim([t(1) t(end)]);
        ylim(yl);
        if m == 3
            xlabel('Time (s)');
        end

        if isempty(hFault)
            legend([hFD, hTh], {'Detection Function', 'Threshold'}, 'Location', 'best');
        else
            legend([hFault, hFD, hTh], {'Fault Period', 'Detection Function', ...
                'Threshold'}, 'Location', 'best');
        end
        grid on;
        hold off;
    end

    saveas(gcf, sprintf('cond%d_detection.fig', cond_num));
    saveas(gcf, sprintf('cond%d_detection.png', cond_num));

    fault_idx = find(fault_mask);
    if ~isempty(fault_idx)
        figure('Name', sprintf('Condition %d - Isolation', cond_num), ...
            'Position', [100 100 1200 500]);

        FIs = {FI_glt, FI_wglt, FI_swglt};
        for m = 1:3
            subplot(1, 3, m);
            FI_mean = mean(FIs{m}(fault_idx, :), 1);
            bar(1:9, FI_mean, 'FaceColor', colors_m{m}, 'EdgeColor', 'k');
            set(gca, 'XTick', 1:9, 'XTickLabel', sensor_labels, 'XTickLabelRotation', 45);
            ylabel('FI (mean during fault)');
            title(sprintf('%s Isolation', methods{m}));
            grid on;

            [~, max_idx] = max(FI_mean);
            hold on;
            bar(max_idx, FI_mean(max_idx), 'FaceColor', [1 0.3 0.3], 'EdgeColor', 'k');
            hold off;
        end

        sgtitle(cond_title);
        saveas(gcf, sprintf('cond%d_isolation.fig', cond_num));
        saveas(gcf, sprintf('cond%d_isolation.png', cond_num));
    end

    if nargin < 16 || isempty(results_swglt)
        return;
    end

    figure('Name', sprintf('Condition %d - SWGLT Diagnostics', cond_num), ...
        'Position', [120 120 1200 650]);

    subplot(2, 1, 1);
    FI_sorted_all = sort(results_swglt.FI, 2, 'descend');
    rho_series = FI_sorted_all(:, 1) ./ max(FI_sorted_all(:, 2), eps);
    rho_series(results_swglt.FD_energy <= results_swglt.T_adaptive) = NaN;
    semilogy(t, rho_series, 'b-', 'LineWidth', 0.8);
    hold on;
    yline(cfg.rho_threshold, 'r--', 'LineWidth', 1.5, 'Label', '\rho_{th}');
    xlabel('Time (s)');
    ylabel('Isolation Margin Ratio \rho');
    title('Isolation Margin Ratio (FI_{(1)} / FI_{(2)})');
    legend('Margin ratio', 'Threshold');
    grid on;

    subplot(2, 1, 2);
    post = results_swglt.posterior;
    [~, top1_idx] = max(post, [], 2);
    plot(t, post(:, true_fault_sensor), 'g-', 'LineWidth', 1.2);
    hold on;

    competitors = top1_idx;
    competitors(competitors == true_fault_sensor) = NaN;
    if any(~isnan(competitors))
        comp_sensor = mode(competitors(~isnan(competitors)));
        plot(t, post(:, comp_sensor), 'r-', 'LineWidth', 0.8);
    else
        comp_sensor = true_fault_sensor;
    end

    yline(cfg.P_isol, 'k--', 'LineWidth', 1.0, 'Label', 'P_{isol}');
    xlabel('Time (s)');
    ylabel('Posterior Probability');
    title('Bayesian Isolation Posterior');
    legend(sprintf('Sensor %d (true)', true_fault_sensor), ...
           sprintf('Sensor %d (competitor)', comp_sensor), 'Threshold');
    grid on;

    saveas(gcf, sprintf('cond%d_swglt_diagnostics.fig', cond_num));
    saveas(gcf, sprintf('cond%d_swglt_diagnostics.png', cond_num));
end
