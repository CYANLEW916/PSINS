function Z_faulty = inject_fault(Z, t, fault_cfg, fault_type)
% INJECT_FAULT  Inject sensor faults into measurement vectors.
%   Z_faulty = INJECT_FAULT(Z, t, fault_cfg, fault_type) adds fault signals
%   to the specified sensor channel in the N x 9 measurement matrix Z.
%
%   Inputs:
%     Z          - N x 9 measurement matrix (clean + noise)
%     t          - N x 1 time vector (s)
%     fault_cfg  - fault configuration structure from config.m
%     fault_type - 'hard' for step faults, 'soft' for ramp faults
%
%   Output:
%     Z_faulty   - N x 9 measurement matrix with injected faults
%
% See also  config, main_simulation.

    Z_faulty = Z;
    idx = fault_cfg.sensor_idx;

    if strcmp(fault_type, 'hard')
        % Hard fault: step bias during specified intervals
        for fi = 1:size(fault_cfg.intervals, 1)
            t_start = fault_cfg.intervals(fi, 1);
            t_end   = fault_cfg.intervals(fi, 2);
            mag     = fault_cfg.magnitudes(fi);

            mask = (t >= t_start) & (t <= t_end);
            Z_faulty(mask, idx) = Z_faulty(mask, idx) + mag;
        end

    elseif strcmp(fault_type, 'soft')
        % Soft fault: linear ramp during specified interval
        t_start = fault_cfg.interval(1);
        t_end   = fault_cfg.interval(2);
        rate    = fault_cfg.rate;
        t_ref   = fault_cfg.t_start_ref;

        mask = (t >= t_start) & (t <= t_end);
        t_fault = t(mask);
        ramp = rate * (t_fault - t_ref);
        Z_faulty(mask, idx) = Z_faulty(mask, idx) + ramp;
    end

end
