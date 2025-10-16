function [err, lat, errStats] = avpcmpplot(avp0, varargin)
% AVPs comparison & errors plot.
%
% Prototype: err = avpcmpplot(avp0, varargin)
% Inputs: avp0 - AVP reference, may be [att], [att,vn], [att,vn,pos]
%                or [vn,pos]
%         varargin - the last input parameter may be comparison type or by
%                    default value
% Output: err - AVP error
%         lat - latitude
%         errStats - statistics (mean & RMS) for available error components
%          
% See also  avpcmp, inserrplot, posflucmpplot, kfplot, gpsplot, imuplot, avperrstd, errstplot, pos2dplot.

% Copyright(c) 2009-2024, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 06/10/2013, 11/10/2024
global glv
    [m, n] = size(avp0);
    if n==1, n=m; m=1; end
    if m==1  % static base avp0, row record
        avp0 = [repmat(avp0', length(varargin{1}), 1), varargin{1}(:,end)];
        n = n+1;
    end
    phi_mu = 'phi';
    errStats = struct();
    lat = [];
    if ischar(varargin{end})
        if strcmp(varargin{end},'phi')==1, phi_mu = 'phi'; varargin = varargin(1:end-1);
        elseif strcmp(varargin{end},'mu')==1, phi_mu = 'mu'; varargin = varargin(1:end-1);
        elseif strcmp(varargin{end},'datt')==1, phi_mu = 'datt'; varargin = varargin(1:end-1); end
    end
    if ischar(varargin{end}),  ptype = varargin{end};  varargin = varargin(1:end-1);
    else
        if n<3,  	ptype = 'y';
        elseif n<6,	ptype = 'a';
        elseif n<9,	
            if max(abs(avp0(:,4)))>pi/2, ptype = 'av'; else, ptype = 'vp'; end
        else      	ptype = 'avp';
        end
    end
    %%
    plotReferenceFigures = false;
    if plotReferenceFigures
        myfig;
%     if size(avp0,2)>10, avp0=avp0(:,[1:9,end]); end
        switch ptype
            case {'A', 'a', 'y'},
                avp = avp0; t = avp(:,end); if size(avp,2)<3, avp0=[zeros(length(avp),2),avp]; avp=avp0; end
                subplot(221), plot(t, avp(:,1:2)/glv.deg), xygo('pr');  legend('Pitch Ref.', 'Roll Ref.');
                subplot(223), plot(t, avp(:,3)/glv.deg), xygo('y');
            case 'av',
                avp = avp0; t = avp(:,end);
                subplot(221), plot(t, avp(:,1:3)/glv.deg), xygo('att');  legend('Pitch Ref.', 'Roll Ref.', 'Yaw Ref.');
                subplot(223), plot(t, [avp(:,4:6),sqrt(avp(:,4).^2+avp(:,5).^2+avp(:,6).^2)]); xygo('V');
            case 'avp',
                avp = avp0; t = avp(:,end);
                subplot(321), plot(t, avp(:,1:3)/glv.deg), xygo('att');  legend('Pitch Ref.', 'Roll Ref.', 'Yaw Ref.');
                subplot(323), plot(t, [avp(:,4:6),sqrt(avp(:,4).^2+avp(:,5).^2+avp(:,6).^2)]); xygo('V');
                subplot(325), plot(t, [[avp(:,7)-avp0(1,7),(avp(:,8)-avp0(1,8))*cos(avp0(1,7))]*glv.Re,avp(:,9)-avp0(1,9)]); xygo('DP');
            case 'vp',
                avp = avp0; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                if size(avp0,2)<10, avp0 = [zeros(length(avp0),3), avp0]; end
    %             subplot(221), plot(t, [avp(:,4:6),sqrt(avp(:,4).^2+avp(:,5).^2+avp(:,6).^2)]); xygo('V');
                subplot(221), plot(t, avp(:,4:6)); xygo('V');   legend('V_E Ref.', 'V_N Ref.', 'V_U Ref.');
                subplot(223), plot(t, [[avp(:,7)-avp0(1,7),(avp(:,8)-avp0(1,8))*cos(avp0(1,7))]*glv.Re,avp(:,9)-avp0(1,9)]); xygo('DP');
            case 'v',
                avp = avp0; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                if size(avp0,2)<10, avp0 = [zeros(length(avp0),3), avp0]; end
    %             subplot(121), plot(t, [avp(:,4:6),sqrt(avp(:,4).^2+avp(:,5).^2+avp(:,6).^2)]); xygo('V');
                subplot(121), plot(t, avp(:,4:6)); xygo('V');
            case 've',
                avp = avp0; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                if size(avp0,2)<10, avp0 = [zeros(length(avp0),3), avp0]; end
                subplot(121), plot(t, avp(:,4)); xygo('V');
            case 'vn',
                avp = avp0; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                if size(avp0,2)<10, avp0 = [zeros(length(avp0),3), avp0]; end
                subplot(121), plot(t, avp(:,5)); xygo('V');
            case 'vu',
                avp = avp0; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                if size(avp0,2)<10, avp0 = [zeros(length(avp0),3), avp0]; end
                subplot(121), plot(t, avp(:,6)); xygo('V');
            case 'vb',
                avp = avp0; t = avp(:,end);
                vb0 = amulvBatch(avp(:,1:3), avp(:,[4:6,end]));
                subplot(121), plot(t, vb0(:,1:3)); xygo('v^b / m/s');
            case 'p',
                avp = avp0; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                if size(avp0,2)<10, avp0 = [zeros(length(avp0),3), avp0]; end
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                if size(avp0,2)<10, avp0 = [zeros(length(avp0),3), avp0]; end
                subplot(121), plot(t, [[avp(:,7)-avp0(1,7),(avp(:,8)-avp0(1,8))*cos(avp0(1,7))]*glv.Re,avp(:,9)-avp0(1,9)]); xygo('DP');
            case 'ed',
                avp = avp0; t = avp(:,end);
                subplot(211), plot(t, avp(:,end-6:end-4)/glv.dph);
                subplot(212), plot(t, avp(:,end-3:end-1)/glv.ug);
        end
    else
        if strcmp(ptype, 'vb')
            vb0 = amulvBatch(avp0(:,1:3), avp0(:,[4:6,end]));
        end
    end
    kk = length(varargin);
    str = '-.--: ';
%     for k=1:kk
%         if size(varargin{k},2)>10, varargin{k}=varargin{k}(:,[1:9,end]); end
%     end
    switch ptype
        case 'A',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
%                subplot(121), hold on, plot(t, avp(:,1:2)/glv.deg, strk, 'LineWidth',2), xygo('pr');
                if plotReferenceFigures
                    subplot(221), hold on, plot(t, avp(:,1)/glv.deg, 'r', t, avp(:,2)/glv.deg, 'm'), xygo('pr');
                    subplot(222), hold on, plot(t, avp(:,3)/glv.deg, 'r'), xygo('y');
                end
%                err = avpcmp(avp, avp0, 'mu'); t = err(:,end);
%                 subplot(122), hold on, plot(t, err(:,1:2)/glv.min, strk, 'LineWidth',2); xygo('mu'); mylegend('mux','muy');
            end
        case {'a', 'y'},
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);  if size(avp,2)<3, avp=[zeros(length(avp),2),avp]; end
                if plotReferenceFigures
                    subplot(221), hold on, plot(t, avp(:,1:2)/glv.deg, strk, 'LineWidth',2), xygo('pr');
                    subplot(223), hold on, plot(t, avp(:,3)/glv.deg, strk, 'LineWidth',2), xygo('y');
                end
                err = avpcmp(avp(:,[1:3,end]), avp0(:,[1:3,end]), phi_mu); t = err(:,end);
                attLabels = attitudeErrorLabels(phi_mu);
                for idx = 1:3
                    plotErrorComponent(ptype, attLabels{idx}, t, err(:,idx)/glv.min, strk, 'Attitude Error', k, kk, groupedErrorFigure('attitude', idx));
                    recordStats('attitude', attLabels{idx}, err(:,idx)/glv.min, k, 'arcmin');
                end
            end
        case 'av',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
                if plotReferenceFigures
                    subplot(221), hold on, plot(t, avp(:,1:3)/glv.deg, strk, 'LineWidth',2), xygo('att');
                    subplot(223), hold on, plot(t, [avp(:,4:6),sqrt(avp(:,4).^2+avp(:,5).^2+avp(:,6).^2)], strk, 'LineWidth',2); xygo('V');
                end
                err = avpcmp(avp(:,[1:6,end]), avp0(:,[1:6,end]), phi_mu); t = err(:,end);
                attLabels = attitudeErrorLabels(phi_mu);
                for idx = 1:3
                    plotErrorComponent(ptype, attLabels{idx}, t, err(:,idx)/glv.min, strk, 'Attitude Error', k, kk, groupedErrorFigure('attitude', idx));
                    recordStats('attitude', attLabels{idx}, err(:,idx)/glv.min, k, 'arcmin');
                end
                velLabels = {'dVE','dVN','dVU'};
                for idx = 1:3
                    plotErrorComponent(ptype, velLabels{idx}, t, err(:,idx+3), strk, 'Velocity Error', k, kk, groupedErrorFigure('velocity', idx));
                    recordStats('velocity', velLabels{idx}, err(:,idx+3), k, 'm/s');
                end
            end
        case 'avp',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
                if plotReferenceFigures
                    subplot(321), hold on, plot(t, avp(:,1:3)/glv.deg, strk, 'LineWidth',2), xygo('att');
                    subplot(323), hold on, plot(t, [avp(:,4:6),sqrt(avp(:,4).^2+avp(:,5).^2+avp(:,6).^2)], strk, 'LineWidth',2); xygo('V');
                    subplot(325), hold on, plot(t, [[avp(:,7)-avp0(1,7),(avp(:,8)-avp0(1,8))*cos(avp0(1,7))]*glv.Re,avp(:,9)-avp0(1,9)], strk, 'LineWidth',2); xygo('DP');
                end
                [err,i1,i0] = avpcmp(avp(:,[1:9,end]), avp0(:,[1:9,end]), phi_mu); t = err(:,end);  lat = avp0(i0,7);
                attLabels = attitudeErrorLabels(phi_mu);
                for idx = 1:3
                    plotErrorComponent(ptype, attLabels{idx}, t, err(:,idx)/glv.min, strk, 'Attitude Error', k, kk, groupedErrorFigure('attitude', idx));
                    recordStats('attitude', attLabels{idx}, err(:,idx)/glv.min, k, 'arcmin');
                end
                velLabels = {'dVE','dVN','dVU'};
                for idx = 1:3
                    plotErrorComponent(ptype, velLabels{idx}, t, err(:,idx+3), strk, 'Velocity Error', k, kk, groupedErrorFigure('velocity', idx));
                    recordStats('velocity', velLabels{idx}, err(:,idx+3), k, 'm/s');
                end
                posData = {err(:,7)*glv.Re, err(:,8).*cos(lat)*glv.Re, err(:,9)};
                posLabels = {'dlat','dlon','dH'};
                for idx = 1:3
                    plotErrorComponent(ptype, posLabels{idx}, t, posData{idx}, strk, 'Position Error', k, kk, groupedErrorFigure('position', idx));
                    recordStats('position', posLabels{idx}, posData{idx}, k, 'm');
                end
            end
        case 'avpdist',  % subplot(326), pos_err vs. distance
            avp = varargin{1};
            close(gcf);
            err = avpcmpplot(avp0, avp, phi_mu);
            pos = avp(:,[7:9,end]);
            [RMh, clRNh] = RMRN(pos);
            dpos = [zeros(1,3);diff(pos(:,1:3))];
            dxyz = [dpos(:,2).*clRNh, dpos(:,1).*RMh, dpos(:,3)];
            dist = cumsum(normv(dxyz));
            dist = interp1(pos(:,end), dist, err(:,end));
            t1 = get(gca,'xlim');
            distFigData = {err(:,7)*glv.Re, err(:,8)*cos(avp(1,7))*glv.Re, err(:,9)};
            posLabels = {'dlat','dlon','dH'};
            for idx = 1:3
                ax = plotErrorComponent('avpdist', posLabels{idx}, dist/1000, distFigData{idx}, 'k-', 'Position Error vs Distance', 1, 1);
                x1 = (t1(1)-err(1,end))/(err(1,end)-err(end,end))*(dist(1,1)-dist(end,1))+dist(1,1);
                x2 = (t1(end)-err(1,end))/(err(1,end)-err(end,end))*(dist(1,1)-dist(end,1))+dist(1,1);
                xlim(ax, [x1,x2]/1000);
                xlabel(ax, labeldef('dist / km'));
                recordStats('positionDist', posLabels{idx}, distFigData{idx}, 1, 'm');
            end
        case 'vp',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
%                 subplot(221), hold on, plot(t, [avp(:,4:6),sqrt(avp(:,4).^2+avp(:,5).^2+avp(:,6).^2)], strk, 'LineWidth',2); xygo('V');
                if plotReferenceFigures
                    subplot(221), hold on, plot(t, avp(:,4:6), strk, 'LineWidth',2); xygo('V');
                    subplot(223), hold on, plot(t, [[avp(:,7)-avp0(1,7),(avp(:,8)-avp0(1,8))*cos(avp0(1,7))]*glv.Re,avp(:,9)-avp0(1,9)], strk, 'LineWidth',2); xygo('DP');
                end
                [err,i1,i0] = avpcmp(avp, avp0, 'noatt'); t = err(:,end);  lat = avp0(i0,7);
                velLabels = {'dVE','dVN','dVU'};
                for idx = 1:3
                    plotErrorComponent(ptype, velLabels{idx}, t, err(:,idx+3), strk, 'Velocity Error', k, kk, groupedErrorFigure('velocity', idx));
                    recordStats('velocity', velLabels{idx}, err(:,idx+3), k, 'm/s');
                end
                posData = {err(:,7)*glv.Re, err(:,8).*cos(lat)*glv.Re, err(:,9)};
                posLabels = {'dlat','dlon','dH'};
                for idx = 1:3
                    plotErrorComponent(ptype, posLabels{idx}, t, posData{idx}, strk, 'Position Error', k, kk, groupedErrorFigure('position', idx));
                    recordStats('position', posLabels{idx}, posData{idx}, k, 'm');
                end
            end
        case 'v',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
%                 subplot(121), hold on, plot(t, [avp(:,4:6),sqrt(avp(:,4).^2+avp(:,5).^2+avp(:,6).^2)], strk, 'LineWidth',2); xygo('V');
                if plotReferenceFigures
                    subplot(121), hold on, plot(t, avp(:,4:6), strk, 'LineWidth',2); xygo('V');
                end
                err = avpcmp(avp, avp0, 'noatt'); t = err(:,end);
                velLabels = {'dVE','dVN','dVU'};
                for idx = 1:3
                    plotErrorComponent(ptype, velLabels{idx}, t, err(:,idx+3), strk, 'Velocity Error', k, kk, groupedErrorFigure('velocity', idx));
                    recordStats('velocity', velLabels{idx}, err(:,idx+3), k, 'm/s');
                end
            end
        case 've',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                err = avpcmp(avp, avp0, 'noatt'); t = err(:,end);
                plotErrorComponent(ptype, 'dVE', t, err(:,4), strk, 'Velocity Error', k, kk, groupedErrorFigure('velocity', 1));
                recordStats('velocity', 'dVE', err(:,4), k, 'm/s');
            end
        case 'vn',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                err = avpcmp(avp, avp0, 'noatt'); t = err(:,end);
                plotErrorComponent(ptype, 'dVN', t, err(:,5), strk, 'Velocity Error', k, kk, groupedErrorFigure('velocity', 2));
                recordStats('velocity', 'dVN', err(:,5), k, 'm/s');
            end
        case 'vu',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                err = avpcmp(avp, avp0, 'noatt'); t = err(:,end);
                plotErrorComponent(ptype, 'dVU', t, err(:,6), strk, 'Velocity Error', k, kk, groupedErrorFigure('velocity', 3));
                recordStats('velocity', 'dVU', err(:,6), k, 'm/s');
            end
        case 'vb',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
                vb = amulvBatch(avp(:,1:3),avp(:,[4:6,end]));
                err = avpcmp(vb, vb0, 'noatt'); t = err(:,end);
                vbLabels = {'dVR','dVF','dVUp'};
                for idx = 1:3
                    plotErrorComponent(ptype, vbLabels{idx}, t, err(:,idx), strk, 'Body Velocity Error', k, kk, groupedErrorFigure('velocity', idx));
                    recordStats('bodyVelocity', vbLabels{idx}, err(:,idx), k, 'm/s');
                end
            end
        case 'p',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                if size(avp,2)<10, avp = [zeros(length(avp),3), avp]; end
                [err,i1,i0] = avpcmp(avp, avp0, 'noatt'); t = err(:,end);  lat = avp0(i0,7);
                posData = {err(:,7)*glv.Re, err(:,8).*cos(lat)*glv.Re, err(:,9)};
                posLabels = {'dlat','dlon','dH'};
                for idx = 1:3
                    plotErrorComponent(ptype, posLabels{idx}, t, posData{idx}, strk, 'Position Error', k, kk, groupedErrorFigure('position', idx));
                    recordStats('position', posLabels{idx}, posData{idx}, k, 'm');
                end
            end
        case 'ed',
            for k=1:kk
                strk = str(k*2-1:k*2);
                avp = varargin{k}; t = avp(:,end);
                if plotReferenceFigures
                    subplot(211), hold on, plot(t, avp(:,end-6:end-4)/glv.dph, strk, 'LineWidth',2); xygo('eb');
                    subplot(212), hold on, plot(t, avp(:,end-3:end-1)/glv.ug, strk, 'LineWidth',2); xygo('db');
                end
            end
    end
    
    function ax = plotErrorComponent(ptypeStr, compLabel, xData, yData, styleSpec, titleStr, seriesIdx, totalSeries, groupInfo)
        if nargin < 5 || isempty(styleSpec)
            styleSpec = '-';
        end
        if nargin < 6
            titleStr = '';
        end
        if nargin < 7
            seriesIdx = [];
        end
        if nargin < 8
            totalSeries = [];
        end
        if nargin < 9 || isempty(groupInfo)
            groupInfo = struct();
        end
        if ~isempty(groupInfo)
            tag = groupInfo.figureTag;
        else
            tag = ['avpcmpplot_', ptypeStr, '_', compLabel];
        end
        [ax, axCreated] = ensureErrorAxis(tag, compLabel, titleStr, groupInfo);
        shouldReset = axCreated || (~isempty(seriesIdx) && seriesIdx == 1);
        if shouldReset
            cla(ax);
            xygo(compLabel);
            if ~isempty(titleStr)
                title(ax, labeldef(titleStr));
            end
        end
        hold(ax, 'on');
        grid(ax, 'on');
        if isempty(seriesIdx) || isempty(totalSeries) || totalSeries <= 1
            plot(ax, xData, yData, styleSpec, 'LineWidth',2);
            if ~isempty(totalSeries) && totalSeries <= 1
                legend(ax, 'off');
            end
        else
            plot(ax, xData, yData, styleSpec, 'LineWidth',2, ...
                'DisplayName', sprintf('Run %d', seriesIdx));
            legend(ax, 'show', 'Location', 'best');
        end
    end

    function [ax, created] = ensureErrorAxis(tag, labelStr, titleStr, groupInfo)
        if nargin < 4
            groupInfo = struct();
        end
        created = false;
        if isempty(groupInfo)
            figTag = tag;
            axesKey = 'mainAxes';
            figName = [labelStr, ' Error'];
        else
            figTag = tag;
            axesKey = sprintf('axes_%d', groupInfo.subplotIndex);
            if isfield(groupInfo, 'figureName') && ~isempty(groupInfo.figureName)
                figName = groupInfo.figureName;
            else
                figName = [titleStr, ' Components'];
            end
        end
        fig = findobj(0, 'Type', 'figure', 'Tag', figTag);
        if isempty(fig) || ~ishandle(fig)
            fig = figure('Name', figName, 'Tag', figTag);
            created = true;
        else
            figure(fig);
        end
        if isempty(groupInfo)
            ax = getappdata(fig, axesKey);
            if isempty(ax) || ~ishandle(ax)
                clf(fig);
                ax = axes('Parent', fig);
                setappdata(fig, axesKey, ax);
                created = true;
            end
        else
            ax = getappdata(fig, axesKey);
            if isempty(ax) || ~ishandle(ax)
                if ~isfield(groupInfo, 'numRows') || isempty(groupInfo.numRows)
                    groupInfo.numRows = 3;
                end
                if ~isfield(groupInfo, 'numCols') || isempty(groupInfo.numCols)
                    groupInfo.numCols = 1;
                end
                ax = subplot(groupInfo.numRows, groupInfo.numCols, groupInfo.subplotIndex, 'Parent', fig);
                setappdata(fig, axesKey, ax);
                created = true;
            end
        end
        axes(ax);
    end

    function info = groupedErrorFigure(groupName, subplotIndex)
        if nargin < 2
            subplotIndex = 1;
        end
        switch lower(groupName)
            case 'attitude'
                info.figureTag = 'avpcmpplot_group_attitude';
                info.figureName = 'Attitude Error Components';
            case {'velocity','bodyvelocity'}
                info.figureTag = 'avpcmpplot_group_velocity';
                info.figureName = 'Velocity Error Components';
            case 'position'
                info.figureTag = 'avpcmpplot_group_position';
                info.figureName = 'Position Error Components';
            otherwise
                info.figureTag = ['avpcmpplot_group_', groupName];
                info.figureName = [capitalizeFirst(groupName), ' Error Components'];
        end
        info.numRows = 3;
        info.numCols = 1;
        info.subplotIndex = subplotIndex;
    end

    function labels = attitudeErrorLabels(mode)
        if isempty(mode)
            mode = 'phi';
        end
        switch mode(1)
            case 'm'
                labels = {'mux','muy','muz'};
            case 'd'
                labels = {'dpch','drll','dyaw'};
            otherwise
                labels = {'phiE','phiN','phiU'};
        end
    end

    function recordStats(group, component, data, runIdx, unitStr)
        if nargin < 5
            unitStr = '';
        end
        if ~isfield(errStats, group)
            errStats.(group) = struct([]);
        end
        data = data(:);
        data = data(~isnan(data));
        if isempty(data)
            meanVal = NaN;
            rmsVal = NaN;
        else
            meanVal = mean(data);
            rmsVal = sqrt(mean(data.^2));
        end
        statsEntry = struct('mean', meanVal, 'rms', rmsVal, 'unit', unitStr);
        if numel(errStats.(group)) < runIdx || isempty(errStats.(group))
            errStats.(group)(runIdx) = struct();
        end
        errStats.(group)(runIdx).(component) = statsEntry;
        if isempty(unitStr)
            unitFmt = '';
        else
            unitFmt = [' ', unitStr];
        end
        fprintf('%s error (%s) run %d: mean = %.6g%s, RMS = %.6g%s\n', ...
            capitalizeFirst(group), component, runIdx, meanVal, unitFmt, rmsVal, unitFmt);
    end

    function out = capitalizeFirst(str)
        if isempty(str)
            out = str;
            return;
        end
        out = [upper(str(1)), str(2:end)];
    end

end
