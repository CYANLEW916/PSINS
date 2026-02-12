function stats = evaluate_performance(FD, threshold, fault_mask, method_name, is_soft, t)
% EVALUATE_PERFORMANCE  Compute FDI statistical performance metrics.
%   stats = EVALUATE_PERFORMANCE(FD, threshold, fault_mask, method_name, is_soft, t)
%   computes fault detection rate, false alarm rate, accuracy, and detection delay.
%
%   Inputs:
%     FD          - N x 1 detection function values
%     threshold   - scalar or N x 1 detection threshold
%     fault_mask  - N x 1 logical vector (true during fault periods)
%     method_name - string name of the method for display
%     is_soft     - true for soft faults (compute detection delay)
%     t           - N x 1 time vector (s)
%
%   Output:
%     stats       - structure with fields: FDR, FAR, Accuracy, delay
%
% See also  main_simulation, fdi_glt, fdi_wglt, fdi_swglt.

    N = length(FD);

    % Detection decisions
    if isscalar(threshold)
        detected = FD > threshold;
    else
        detected = FD > threshold;
    end

    % Fault period indices
    fault_idx   = fault_mask;
    normal_idx  = ~fault_mask;

    n_fault  = sum(fault_idx);
    n_normal = sum(normal_idx);

    % True positives: correctly detected faults
    TP = sum(detected & fault_idx);

    % False positives: false alarms during normal operation
    FP = sum(detected & normal_idx);

    % True negatives: correctly identified normal periods
    TN = sum(~detected & normal_idx);

    % False negatives: missed faults
    FN = sum(~detected & fault_idx);

    % Fault Detection Rate (FDR)
    if n_fault > 0
        stats.FDR = TP / n_fault;
    else
        stats.FDR = 0;
    end

    % False Alarm Rate (FAR)
    if n_normal > 0
        stats.FAR = FP / n_normal;
    else
        stats.FAR = 0;
    end

    % Accuracy
    stats.Accuracy = (TP + TN) / N;

    % Detection delay (only for soft faults)
    stats.delay = NaN;
    if is_soft && n_fault > 0
        % Find the start of the fault period
        fault_start_idx = find(fault_idx, 1, 'first');
        t_fault_start = t(fault_start_idx);

        % Find the first detection within the fault period
        detected_in_fault = detected & fault_idx;
        first_detection_idx = find(detected_in_fault, 1, 'first');

        if ~isempty(first_detection_idx)
            stats.delay = t(first_detection_idx) - t_fault_start;
        else
            stats.delay = Inf;  % never detected
        end
    end

    fprintf('  %s: FDR=%.2f%%, FAR=%.2f%%, Accuracy=%.2f%%', ...
        method_name, stats.FDR*100, stats.FAR*100, stats.Accuracy*100);
    if is_soft
        if isinf(stats.delay)
            fprintf(', Delay=N/A (not detected)');
        else
            fprintf(', Delay=%.1fs', stats.delay);
        end
    end
    fprintf('\n');

end
