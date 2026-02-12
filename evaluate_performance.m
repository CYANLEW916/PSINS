function stats = evaluate_performance(FD, threshold, fault_mask, method_name, is_soft, t)
% EVALUATE_PERFORMANCE  Compute FDI statistical performance metrics.
%   stats = EVALUATE_PERFORMANCE(FD, threshold, fault_mask, method_name, is_soft, t)
%
%   Metrics:
%     FDR      = correctly detected fault samples / total fault samples
%     FAR      = false alarms in normal period  / total normal samples
%     Accuracy = (TP + TN) / total samples
%     Delay    = time from fault onset to first detection (soft faults only)
%
%   Inputs:
%     FD          - N x 1 detection function values
%     threshold   - scalar detection threshold
%     fault_mask  - N x 1 logical vector (true = fault period)
%     method_name - string label for display
%     is_soft     - true for soft faults (compute detection delay)
%     t           - N x 1 time vector (s)
%
%   Output:
%     stats       - structure with fields FDR, FAR, Accuracy, delay
%
% See also  main_simulation, fdi_glt, fdi_wglt, fdi_swglt.

    N = length(FD);
    detected = FD > threshold;

    n_fault  = sum(fault_mask);
    n_normal = sum(~fault_mask);

    TP = sum(detected & fault_mask);     % correct detections
    FP = sum(detected & ~fault_mask);    % false alarms
    TN = sum(~detected & ~fault_mask);   % correct normals

    % FDR
    if n_fault > 0
        stats.FDR = TP / n_fault;
    else
        stats.FDR = 0;
    end

    % FAR
    if n_normal > 0
        stats.FAR = FP / n_normal;
    else
        stats.FAR = 0;
    end

    % Accuracy
    stats.Accuracy = (TP + TN) / N;

    % Detection delay (soft faults only)
    stats.delay = NaN;
    if is_soft && n_fault > 0
        fault_start_idx = find(fault_mask, 1, 'first');
        t_fault_start = t(fault_start_idx);

        first_det = find(detected & fault_mask, 1, 'first');
        if ~isempty(first_det)
            stats.delay = t(first_det) - t_fault_start;
        else
            stats.delay = Inf;
        end
    end

    % Console output
    fprintf('  %s: FDR=%.2f%%, FAR=%.2f%%, Accuracy=%.2f%%', ...
        method_name, stats.FDR*100, stats.FAR*100, stats.Accuracy*100);
    if is_soft
        if isinf(stats.delay)
            fprintf(', Delay=N/A');
        else
            fprintf(', Delay=%.1fs', stats.delay);
        end
    end
    fprintf('\n');

end
