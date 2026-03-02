function [Z, meta] = fdi_generate_redundant_measurements(ref3, sigmaIns, sigmaIsis, faultCfg)
%FDI_GENERATE_REDUNDANT_MEASUREMENTS Create 2 INS + 1 ISIS heterogeneous measurements.
% Inputs:
%   ref3      - Reference signal (N x 3), e.g., gyro rate or specific force.
%   sigmaIns  - INS per-axis 1-sigma noise standard deviation.
%   sigmaIsis - ISIS per-axis 1-sigma noise standard deviation.
%   faultCfg  - Optional fault configuration struct. Supported modes:
%               (A) single fault with fields:
%                   .sensorIdx, .startIdx, .type('step'|'ramp'), .amplitude
%               (B) campaign/event mode with field:
%                   .events: struct array with fields
%                   .sensorIdx, .startIdx, .endIdx, .type, .amplitude
% Output:
%   Z         - Redundant measurement matrix (N x 9).
%   meta      - Struct including noise sigma vector and injected fault profile.

    if nargin < 4
        faultCfg = struct();
    end

    N = size(ref3, 1);
    Ztrue = [ref3, ref3, ref3];

    sigmaVec = [repmat(sigmaIns, 6, 1); repmat(sigmaIsis, 3, 1)];
    noise = randn(N, 9) .* sigmaVec';
    Z = Ztrue + noise;

    fault = zeros(N, 9);
    if isfield(faultCfg, 'events')
        events = faultCfg.events;
        for k = 1:numel(events)
            fault = apply_single_event(fault, events(k), N);
        end
    elseif isfield(faultCfg, 'sensorIdx') && isfield(faultCfg, 'startIdx')
        evt.sensorIdx = faultCfg.sensorIdx;
        evt.startIdx = faultCfg.startIdx;
        evt.endIdx = N;
        evt.type = faultCfg.type;
        evt.amplitude = faultCfg.amplitude;
        fault = apply_single_event(fault, evt, N);
    end

    Z = Z + fault;

    meta.sigmaVec = sigmaVec;
    meta.fault = fault;
end

function fault = apply_single_event(fault, evt, N)
%APPLY_SINGLE_EVENT Apply one step/ramp event in [startIdx, endIdx] to one channel.

    idx = evt.sensorIdx;
    k0 = max(1, evt.startIdx);
    if isfield(evt, 'endIdx')
        k1 = min(N, evt.endIdx);
    else
        k1 = N;
    end

    if k0 > k1
        return;
    end

    amp = evt.amplitude;
    ftype = evt.type;
    span = (k0:k1)';

    if strcmpi(ftype, 'step')
        fault(span, idx) = fault(span, idx) + amp;
    elseif strcmpi(ftype, 'ramp')
        fault(span, idx) = fault(span, idx) + amp * (span - k0);
    else
        error('Unknown fault type: %s', ftype);
    end
end
