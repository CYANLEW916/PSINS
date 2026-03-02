function stats = evaluate_performance(FD, threshold, fault_mask, method_name, is_soft, t, varargin)
% EVALUATE_PERFORMANCE  Compute FDI performance metrics including DNI/CIR/MIR.
%   stats = EVALUATE_PERFORMANCE(FD, threshold, fault_mask, method_name,
%   is_soft, t, isolated, iso_status, true_fault_sensor) computes detection,
%   false alarm, accuracy, delay, and isolation robustness metrics.
%
% See also  main_simulation, fdi_swglt.

    N = length(FD);
    detected = FD > threshold;

    n_fault = sum(fault_mask);
    n_normal = sum(~fault_mask);

    TP = sum(detected & fault_mask);
    FP = sum(detected & ~fault_mask);
    TN = sum(~detected & ~fault_mask);

    if n_fault > 0
        stats.FDR = TP / n_fault;
    else
        stats.FDR = 0;
    end

    if n_normal > 0
        stats.FAR = FP / n_normal;
    else
        stats.FAR = 0;
    end

    stats.Accuracy = (TP + TN) / N;

    stats.delay = NaN;
    if is_soft && n_fault > 0
        fault_start_idx = find(fault_mask, 1, 'first');
        first_det = find(detected & fault_mask, 1, 'first');
        if ~isempty(first_det)
            stats.delay = t(first_det) - t(fault_start_idx);
        else
            stats.delay = Inf;
        end
    end

    stats.DNI_rate = NaN;
    stats.CIR = NaN;
    stats.MIR = NaN;

    if numel(varargin) >= 3
        isolated = varargin{1};
        iso_status = varargin{2};
        true_fault_sensor = varargin{3};
        if isempty(isolated) || isempty(iso_status)
            isolated = [];
        end
    else
        isolated = [];
    end

    if ~isempty(isolated)
        fault_present = fault_mask;
        fault_detected_mask = (iso_status >= 1);
        dni_mask = (isolated == -1) & fault_present;
        n_dni = sum(dni_mask);
        denom_dni = max(sum(fault_present & fault_detected_mask), 1);
        stats.DNI_rate = n_dni / denom_dni;

        confirmed_mask = (iso_status == 2) & fault_present;
        correct_mask = confirmed_mask & (isolated == true_fault_sensor);
        stats.CIR = sum(correct_mask) / max(sum(confirmed_mask), 1);
        stats.MIR = 1 - stats.CIR;
    end

    fprintf('  %s: FDR=%.2f%%, FAR=%.2f%%, Accuracy=%.2f%%', ...
        method_name, stats.FDR * 100, stats.FAR * 100, stats.Accuracy * 100);

    if is_soft
        if isinf(stats.delay)
            fprintf(', Delay=N/A');
        else
            fprintf(', Delay=%.1fs', stats.delay);
        end
    end

    if ~isnan(stats.DNI_rate)
        fprintf(', DNI=%.2f%%, CIR=%.2f%%, MIR=%.2f%%', ...
            stats.DNI_rate * 100, stats.CIR * 100, stats.MIR * 100);
    end
    fprintf('\n');
end
