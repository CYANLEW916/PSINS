% TEST_SINS_FDI_SWGLT_TUNING  Example grid sweep for SWGLT isolation tuning.
%   This demo runs main_simulation with several robust-isolation parameter
%   combinations and prints the tested tuples for reproducible analysis.
%
%   Inputs:  none (edits local cfg values in-loop)
%   Outputs: console logs and normal figure outputs from main_simulation flow

clear; clc;

rho_list = [1.1, 1.3, 1.5, 2.0];
post_list = [0.55, 0.65, 0.75, 0.90];
dwell_list = [1, 2, 3];

fprintf('SWGLT tuning grid size: %d\n', numel(rho_list) * numel(post_list) * numel(dwell_list));

for ir = 1:numel(rho_list)
    for ip = 1:numel(post_list)
        for id = 1:numel(dwell_list)
            cfg = config();
            cfg.rho_threshold = rho_list(ir);
            cfg.P_isol = post_list(ip);
            cfg.N_dwell = dwell_list(id);
            cfg.min_isolation_votes = 2;

            fprintf('Test tuple: rho=%.2f, P_isol=%.2f, N_dwell=%d\n', ...
                cfg.rho_threshold, cfg.P_isol, cfg.N_dwell);

            % Run in base workspace to avoid main_simulation 'clear'
            % wiping this script's loop indices (ir/ip/id).
            assignin('base', 'cfg_override', cfg);
            evalin('base', 'main_simulation');
            evalin('base', 'clear cfg_override');
        end
    end
end
