

x = 1:1100;
y1 = sin(x(1:1000)/100);
y2 = sin(x(101:end)/100);

% figure
% plot(y1)
% hold
% plot(y2)


%%
app = InteractiveUI('Data aligning tool',@(app) (AlignUpdate(app)), .2);


app.Data.y1 = y1;
app.Data.x1 = 1:length(y1);
app.Data.y2 = y2;
app.Data.x2 = 1:length(y2);


app.AddSlider('Latency course', 0, [-100 100]) % percent of total time
app.AddSlider('Latency fine', 0, [-100 100]) % 100 samples
app.AddSlider('Latency superfine', 0, [-1000 1000]) % actual samples
app.Open();



function AlignUpdate(app)


maxDuration = max(length(app.Data.x1),length(app.Data.x2));

if ( ~isfield(app.Data, "f") || ~isvalid(app.Data.f))
    % If figure does not exist create it with all the plots and
    % handles to them

    app.Data.f = figure;
    app.Data.h = struct();
    app.Data.h.y1 = plot(app.Data.x1, app.Data.y1);
    hold
    app.Data.h.y2 = plot(app.Data.x2, app.Data.y2);

    set(gca,'xlim',[-maxDuration, maxDuration*2])

    yl = get(gca,'ylim');
    xl = get(gca,'xlim');

    app.Data.htext = text(xl(1),yl(2), sprintf('IDX latency = %d',0),"HorizontalAlignment",'left', 'VerticalAlignment','top');
end

latencyCourse = app.Values.LatencyCourse;
latencyFine = app.Values.LatencyFine;
latencySuperFine = app.Values.LatencySuperfine;

app.Data.idxLatency = round(maxDuration*latencyCourse/100) + latencyFine*100  + latencySuperFine;

set(app.Data.h.y2,'xdata',app.Data.x2 + app.Data.idxLatency);


    app.Data.htext.String = sprintf('IDX latency = %d',app.Data.idxLatency);

end
