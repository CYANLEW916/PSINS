function test_SINS_FDI_SWGLT_tuning()
% TEST_SINS_FDI_SWGLT_TUNING  Grid sweep for SWGLT isolation tuning.
%   TEST_SINS_FDI_SWGLT_TUNING() runs main_simulation with robust-isolation
%   parameter combinations and prints tested tuples for reproducible analysis.
%
%   Inputs:
%     none
%   Outputs:
%     Console logs and figure outputs from main_simulation.

    clc;

    rho_list = [1.1, 1.3, 1.5, 2.0];
    post_list = [0.55, 0.65, 0.75, 0.90];
    dwell_list = [1, 2, 3];

    n_total = numel(rho_list) * numel(post_list) * numel(dwell_list);
    fprintf('SWGLT tuning grid size: %d\n', n_total);

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

                % Run in base workspace so main_simulation can consume
                % cfg_override without modifying its existing script structure.
                assignin('base', 'cfg_override', cfg);
                evalin('base', 'main_simulation');
                evalin('base', 'clear cfg_override');
            end
        end
    end
end
