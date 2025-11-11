function kf = kfinit0(kf, nts)
% Always called by kfinit and initialize the remaining fields of kf.
%
% See also kfinit, kfupdate, kffeedback, psinstypedef.
    kf.nts = nts;
    [kf.m, kf.n] = size(kf.Hk);
    kf.I = eye(kf.n);
    kf.Kk = zeros(kf.n, kf.m);
    kf.measmask = [];            % measurement mask for no update  20/11/2022
    kf.measstop = zeros(kf.m,1); % measurement stop time
    kf.measlost = zeros(kf.m,1); % measurement lost time
    kf.measlog = 0;              % measurement log flag
    if ~isfield(kf, 'xk'),  kf.xk = zeros(kf.n, 1);  end
    if ~isfield(kf, 'Qk'),  kf.Qk = kf.Qt*kf.nts;  end
    if ~isfield(kf, 'Gammak'),  kf.Gammak = 1; kf.l = kf.n;  else, kf.l=size(kf.Gammak,2);  end
    if ~isfield(kf, 'fading'),  kf.fading = 1;  end
    if ~isfield(kf, 'adaptive'),  kf.adaptive = 0;  end
%     if kf.adaptive==1
        if ~isfield(kf, 'b'),  kf.b = 0.5;  end
        if ~isfield(kf, 'beta'),  kf.beta = 1;  end
        if ~isfield(kf, 'Rmin'),  kf.Rmin = 0.01*kf.Rk;  end
        if ~isfield(kf, 'Rmax'),  kf.Rmax = 100*kf.Rk;  end
        if ~isfield(kf, 'Qmin'),  kf.Qmin = 0.01*kf.Qk;  end
        if ~isfield(kf, 'Qmax'),  kf.Qmax = 100*kf.Qk;  end
%     end
    if ~isfield(kf, 'xtau'),  kf.xtau = ones(size(kf.xk))*eps;   end
    if ~isfield(kf, 'T_fb'),  kf.T_fb = 1;   end
    if ~isfield(kf, 'fbstr'),  kf.fbstr = 'avped';  end
    if ~isfield(kf, 'xconstrain'),  kf.xconstrain = 0;  end
    if ~isfield(kf, 'pconstrain'),  kf.pconstrain = 0;  end
    kf.Pmax = (diag(kf.Pxk)+1)*1.0e10;
    kf.Pmin = kf.Pmax*0;
    if ~isfield(kf, 'fd') || ~isstruct(kf.fd)
        kf.fd = struct();
    end
    if ~isfield(kf.fd, 'enable'), kf.fd.enable = 0; end
    if ~isfield(kf.fd, 'chi2Threshold'), kf.fd.chi2Threshold = inf; end
    if ~isfield(kf.fd, 'holdTime'), kf.fd.holdTime = 0; end
    if ~isfield(kf.fd, 'isOutlier'), kf.fd.isOutlier = false(kf.m,1); end
    if ~isfield(kf.fd, 'whitenedResidual'), kf.fd.whitenedResidual = zeros(kf.m,1); end
    if ~isfield(kf.fd, 'nis'), kf.fd.nis = 0; end
    if ~isfield(kf.fd, 'inDwell'), kf.fd.inDwell = false(kf.m,1); end
    if ~isfield(kf.fd, 'slidingEnable'), kf.fd.slidingEnable = 0; end
    if ~isfield(kf.fd, 'stateLag'), kf.fd.stateLag = 5; end
    if ~isfield(kf.fd, 'slidingWindow'), kf.fd.slidingWindow = max(5, kf.fd.stateLag); end
    if ~isfield(kf.fd, 'slidingThreshold'), kf.fd.slidingThreshold = inf; end
    if ~isfield(kf.fd, 'slidingReference'), kf.fd.slidingReference = 1; end
    if ~isfield(kf.fd, 'slidingBuffer'), kf.fd.slidingBuffer = zeros(kf.m, kf.fd.slidingWindow); end
    if ~isfield(kf.fd, 'slidingIndex'), kf.fd.slidingIndex = 1; end
    if ~isfield(kf.fd, 'slidingStat'), kf.fd.slidingStat = zeros(kf.m,1); end
    if ~isfield(kf.fd, 'slidingTriggered'), kf.fd.slidingTriggered = false; end
    if ~isfield(kf.fd, 'slidingOutlierIdx'), kf.fd.slidingOutlierIdx = false(kf.m,1); end
    if ~isfield(kf.fd, 'slidingNis'), kf.fd.slidingNis = zeros(kf.m,1); end
    if ~isfield(kf.fd, 'slidingCusum'), kf.fd.slidingCusum = zeros(kf.m,1); end
    if ~isfield(kf.fd, 'slidingCusumThreshold'), kf.fd.slidingCusumThreshold = inf; end
    if ~isfield(kf.fd, 'slidingCusumDecay'), kf.fd.slidingCusumDecay = 0; end
    if ~isfield(kf.fd, 'slidingCusumTriggered'), kf.fd.slidingCusumTriggered = false(kf.m,1); end
    if ~isfield(kf.fd, 'historyCount'), kf.fd.historyCount = 0; end
    if ~isfield(kf.fd, 'stateHistory') || size(kf.fd.stateHistory,1)~=kf.n
        kf.fd.stateHistory = repmat(kf.xk, 1, kf.fd.stateLag+1);
        kf.fd.historyValid = 1;
    elseif ~isfield(kf.fd, 'historyValid')
        kf.fd.historyValid = min(size(kf.fd.stateHistory,2), kf.fd.stateLag+1);
    end
%     kf.Pykk_1 = kf.Hk*kf.Pxk*kf.Hk'+kf.Rk;
    kf.Pykk_1 = kf.Hk*kf.Pxk*kf.Hk'+0;
    kf.xfb = zeros(kf.n, 1);
%     kf.coef_fb = (1.0-exp(-kf.T_fb./kf.xtau));
%     kf.coef_fb = ar1coefs(kf.T_fb, kf.xtau);
    xtau = kf.xtau;
    xtau(kf.xtau<kf.T_fb) = kf.T_fb;  kf.coef_fb = kf.T_fb./xtau;  %2015-2-22
    