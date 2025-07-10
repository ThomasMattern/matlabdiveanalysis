% Main script body
[resfile, respath] = uigetfile( ...
       {'*.csv','comma separated file (*.csv)'; ...
        '*.txt','tab separated file (*.txt)'; ...
        '*.*',  'All Files (*.*)'}, ...
        'Pick a dive analysis file');
res = [respath resfile];
resdata = result_data2(res);

% Combine Date and Time (EvtOnset & EvtEnd) into full datetime
resdata.EvtOnset = datetime(resdata.Date) + timeofday(resdata.EvtOnset);
resdata.EvtEnd   = datetime(resdata.Date) + timeofday(resdata.EvtEnd);

birds = unique(resdata.birdID);
[indx, tf] = listdlg('PromptString', {'Select bird to visualize dive data:'}, ...
    'SelectionMode','single','ListString',birds,'ListSize',[150,120]);

bird = birds(indx);
subdata = resdata(resdata.birdID == bird, :);

folder = fileparts(res);
raw = fullfile(folder, bird, 'dive.txt');

% Load raw dive data
[sday,smonth,syear,shour,sminute,ssecond,~,depth] = textread(raw, ...
    '%d-%d-%d,%d:%d:%d,%f,%f%*[^\n]', 'headerlines', 1);
depth = baselinetracking(depth);
rawdata = [syear smonth sday shour sminute ssecond depth];

% Plot overview
xaxis = datetime(rawdata(:,1), rawdata(:,2), rawdata(:,3), rawdata(:,4), rawdata(:,5), rawdata(:,6));
yaxis = -rawdata(:,7);

clf;
subplot(4,1,1);
zoom off; pan off; datacursormode off; disableDefaultInteractivity(gca);
ov = plot(xaxis, yaxis);
set(ov, 'ButtonDownFcn', @clickOV, 'PickableParts', 'all', 'HitTest', 'on');
title(bird); ylabel('Depth (m)'); hold on;

assignin('base', 'subdata', subdata);
assignin('base', 'rawdata', rawdata);
assignin('base', 'ov', ov);

% Click handler
function clickOV(hObj, event)
    pt = event.IntersectionPoint;
    if pt(3) == 0
        x = num2ruler(pt(1), hObj.Parent.XAxis);
        assignin('base', 'cDate', x);
        drawDetail('ov', x);
    end
end

function drawDetail(type, cDate)
    disp(['Passed time: ', datestr(cDate, 'dd-mmm-yyyy HH:MM:SS')]);

    sd = evalin('base', 'subdata');
    rd = evalin('base', 'rawdata');

    [~, idx] = min(abs(cDate - sd.EvtOnset));
    tn = sd.TripNumber(idx);
    td = sd(sd.TripNumber == tn, :);

    tripstart = td.EvtOnset(1);
    tripend   = td.EvtEnd(end);

    rdtimes = datetime(rd(:,1), rd(:,2), rd(:,3), rd(:,4), rd(:,5), rd(:,6));
    inTrip = rdtimes >= tripstart & rdtimes <= tripend;

    if ~any(inTrip)
        disp('No matching data in rawdata for this trip.');
        return;
    end

    tripDataX = rdtimes(inTrip);
    tripDataY = -rd(inTrip, 7);

    % === DETAIL PLOT ===
    subplot(4,1,[2 3 4]);
    cla;
    detailPlot = plot(tripDataX, -tripDataY, 'b');
    title(['Trip ', num2str(tn), ' (', datestr(tripstart, 'dd-mmm-yyyy HH:MM'), ')']);
    ylabel('Depth (m)');
    xlabel('Time');
    set(gca, 'YDir', 'reverse');
    axDetail = gca;

    % === LOCK ZOOM TO X ONLY ===
    zoom off;
    zoom(axDetail, 'on');
    z = zoom(gcf);
    set(z, 'Motion', 'horizontal');
    set(z, 'Enable', 'on');
    set(z, 'RightClickAction', 'PostContextMenu');

    % === OVERVIEW PLOT HIGHLIGHT ===
    subplot(4,1,1);
    hold on;

    if evalin('base', 'exist(''tripHighlight'', ''var'')')
        try
            delete(evalin('base','tripHighlight'));
        catch
        end
    end

    tripHighlight = plot(tripDataX, tripDataY, 'r', 'LineWidth', 1.2);
    assignin('base', 'tripHighlight', tripHighlight);

    % === ADD ZOOM LISTENER TO DETAIL PLOT ===
    if evalin('base', 'exist(''xZoomListener'', ''var'')')
        try
            delete(evalin('base', 'xZoomListener'));
        catch
        end
    end
    zoomListener = addlistener(axDetail, 'XLim', 'PostSet', @updateOverviewHighlight);
    assignin('base', 'xZoomListener', zoomListener);

    updateOverviewHighlight();
end


function updateOverviewHighlight(~, ~)
    try
        tripHighlight = evalin('base', 'tripHighlight');
        rawdata = evalin('base', 'rawdata');

        rdtimes = datetime(rawdata(:,1), rawdata(:,2), rawdata(:,3), ...
                           rawdata(:,4), rawdata(:,5), rawdata(:,6));
        rdepths = -rawdata(:,7);

        % Get XLim from the detail plot
        axList = findall(gcf, 'Type', 'axes');
        axDetail = axList(1);  % lowest subplot typically created last
        lims = get(axDetail, 'XLim');

        inZoom = rdtimes >= lims(1) & rdtimes <= lims(2);

        subplot(4,1,1);
        set(tripHighlight, 'XData', rdtimes(inZoom), 'YData', rdepths(inZoom));
    catch ME
        disp(['[Warning] updateOverviewHighlight failed: ', ME.message]);
    end
end
