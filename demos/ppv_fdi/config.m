function cfg = config()
% CONFIG  Parameter configuration for PPV-aided SWGLT fault detection.
%   cfg = CONFIG() returns a structure containing all simulation parameters
%   including sensor noise, fault injection, algorithm settings, and plotting.
%
% See also  main_simulation, fdi_glt, fdi_wglt, fdi_swglt.

    %% Simulation parameters
    cfg.ts = 0.02;           % sampling period (s), 50 Hz
    cfg.rng_seed = 42;       % random seed for reproducibility

    %% Sensor noise characteristics (from paper Table 1)
    % Gyroscope bias instability (deg/h -> rad/s)
    dph2rps = pi / 180 / 3600;  % 1 deg/h in rad/s
    cfg.sigma_ins_gyro  = 0.01 * dph2rps;   % INS1 & INS2 gyro noise std (rad/s)
    cfg.sigma_isis_gyro = 0.1  * dph2rps;   % ISIS gyro noise std (rad/s)

    % Accelerometer bias instability (g -> m/s^2)
    g0 = 9.8;  % gravitational acceleration (m/s^2)
    cfg.sigma_ins_acc  = 1e-4 * g0;         % INS1 & INS2 accel noise std (m/s^2)
    cfg.sigma_isis_acc = 5e-3 * g0;         % ISIS accel noise std (m/s^2)

    % Noise variance vectors (9 sensors each: INS1 xyz, INS2 xyz, ISIS xyz)
    cfg.sigma_gyro = [cfg.sigma_ins_gyro * ones(1,3), ...
                      cfg.sigma_ins_gyro * ones(1,3), ...
                      cfg.sigma_isis_gyro * ones(1,3)];  % 1x9

    cfg.sigma_acc  = [cfg.sigma_ins_acc * ones(1,3), ...
                      cfg.sigma_ins_acc * ones(1,3), ...
                      cfg.sigma_isis_acc * ones(1,3)];   % 1x9

    %% Configuration matrix H (9x3)
    cfg.H = repmat(eye(3), 3, 1);
    cfg.n_sensors = 9;      % total number of sensors per type
    cfg.n_dof = 3;          % degrees of freedom (3 axes)
    cfg.n_parity = cfg.n_sensors - cfg.n_dof;  % parity space dimension = 6

    %% Fault injection parameters
    % Condition 1: INS1 Y-axis gyro hard fault (step fault)
    cfg.fault1.sensor_idx = 2;          % INS1 Y-axis (2nd element in Z_gyro)
    cfg.fault1.type = 'gyro';
    cfg.fault1.intervals = [120 135; 440 455; 595 610];  % time intervals (s)
    cfg.fault1.magnitudes = [0.5; 1.0; 2.0] * dph2rps;  % fault magnitudes (rad/s)

    % Condition 2: INS2 X-axis accelerometer hard fault
    cfg.fault2.sensor_idx = 4;          % INS2 X-axis (4th element in Z_acc)
    cfg.fault2.type = 'acc';
    cfg.fault2.intervals = [100 115; 420 435; 575 590];  % time intervals (s)
    cfg.fault2.magnitudes = [0.005; 0.01; 0.02] * g0;   % fault magnitudes (m/s^2)

    % Condition 3: INS1 Z-axis gyro soft fault (ramp fault)
    cfg.fault3.sensor_idx = 3;          % INS1 Z-axis (3rd element in Z_gyro)
    cfg.fault3.type = 'gyro';
    cfg.fault3.interval = [260 310];    % time interval (s)
    cfg.fault3.rate = 0.02 * dph2rps;   % fault rate (rad/s per second)
    cfg.fault3.t_start_ref = 260;       % reference start time for ramp

    %% GLT parameters
    cfg.PFA = 0.01;                      % probability of false alarm
    cfg.T_D = chi2inv(1 - cfg.PFA, cfg.n_parity);  % detection threshold ~16.81

    %% SWGLT parameters
    cfg.window_length = 25;              % PCA window = 0.5s / 0.02s = 25 samples
    cfg.pca_threshold = 0.85;            % cumulative variance ratio for PCA
    cfg.K_tolerance = 0.01;              % tolerance factor for adaptive threshold

    %% IMU error model parameters for PSINS imuadderr
    % INS1 error model (high precision)
    cfg.ins1_err.eb  = 0.01;    % gyro bias (deg/h)
    cfg.ins1_err.db  = 100;     % accel bias (ug) -> 1e-4 g = 100 ug
    cfg.ins1_err.web = 0.001;   % angular random walk (deg/sqrt(h))
    cfg.ins1_err.wdb = 10;      % velocity random walk (ug/sqrt(Hz))

    % INS2 error model (same as INS1)
    cfg.ins2_err = cfg.ins1_err;

    % ISIS error model (low precision)
    cfg.isis_err.eb  = 0.1;     % gyro bias (deg/h)
    cfg.isis_err.db  = 5000;    % accel bias (ug) -> 5e-3 g = 5000 ug
    cfg.isis_err.web = 0.01;    % angular random walk (deg/sqrt(h))
    cfg.isis_err.wdb = 100;     % velocity random walk (ug/sqrt(Hz))

end
