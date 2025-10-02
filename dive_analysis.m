%INITIALISE RESULT PARAMETERS
TripNumber=  [];
EventNo=     [];
Date=        [];
SurfaceTime= [];
NoFixes=     [];
FirstLat=    [];
FirstLon=    [];
LastLat=     [];
LastLon=     [];
SurfaceDist= [];
SurfaceVelo= [];
EvtOnset=    [];
EvtEnd=      [];
DiveTime=    [];
DiveDistance=[];
EvtMaxDepth= [];
DescDur=     [];
DescVelo=    [];
StartBot=    [];
NoWiggles=   [];
BottomTime=  [];
BottomEnd=   [];
AscDur=      [];
AscVelo=     [];
MeanWigAmp=  [];
MinWigAmp=   [];
MaxWigAmp=   [];
MeanWigVelo= [];
MinWigVelo=  [];
MaxWigVelo=  [];
BenthicDive= [];
DiveSpeed =  [];
latitude =  [];
longitude =  [];
ttf =  [];

% set the analysis script variables from the GPS & SENSOR matrices

% SENSOR = [day month year hour minute second temperature depth];
day = SENSOR(:,1);
month = SENSOR(:,2);
year = SENSOR(:,3);
hour = SENSOR(:,4);
minute = SENSOR(:,5);
second = SENSOR(:,6);
temperature = SENSOR(:,7);
depth = SENSOR(:,8);
sLat = SENSOR(:,9);
sLon = SENSOR(:,10);

% DETERMINE START DATE & TIME
yearcorr=size(num2str(year(1,:)));
if yearcorr(2)<4
    start = datenum((year(1,:)+2000),month(1,:),day(1,:),hour(1,:),minute(1,:),second(1,:));
    timestamp = datenum((year(:,1)+2000),month(:,1),day(:,1),hour(:,1),minute(:,1),second(:,1));
else
    start = datenum((year(1,:)),month(1,:),day(1,:),hour(1,:),minute(1,:),second(1,:));
    timestamp = datenum((year(:,1)),month(:,1),day(:,1),hour(:,1),minute(:,1),second(:,1));
end
yearcorr=[];

% GPS = [day month year hour minute second lat lon hdop speed];

try
    dayGPS = GPS(:,1);
    monthGPS = GPS(:,2);
    yearGPS = GPS(:,3);
    hourGPS = GPS(:,4);
    minuteGPS = GPS(:,5);
    secondGPS = GPS(:,6);
    latGPS = GPS(:,7);
    lonGPS = GPS(:,8);
    
    startGPS = datenum(yearGPS(1,:),monthGPS(1,:),dayGPS(1,:),hourGPS(1,:),minuteGPS(1,:),secondGPS(1,:));
    timestampGPS = datenum(yearGPS(:,1),monthGPS(:,1),dayGPS(:,1),hourGPS(:,1),minuteGPS(:,1),secondGPS(:,1));
    
    % Check if GPS data starts before dive data, and if so, omit from analysis
    out=find(timestampGPS<timestamp(1));
    timestampGPS(out)=[];
    latGPS(out)=[];
    lonGPS(out)=[];
catch
    disp(['no valid GPS data found']);
end




% DETERMINE SAMPLING INTERVAL
if(second(2,:)==0) %account for first line second = 59
    interval = 60-second(1,:);
else
    interval = second(2,:)-second(1,:);
end

if interval<10;
    dateint = 1/100000;
else
    dateint = 1/10000;
end


tit = ['processing...'];
max_count=100;
% BASELINETRACKING
%
% This feature might be used to correct depth values in case the surface
% depth deviates from zero (e.g. surface and on land values are lower than
% zero metres). Under the assumption that this error is constant (e.g. 
% deviations in atmospheric pressure)and thus also applies to positive 
% depth values, the dive values can be corrected by summing the positive 
% of the zero deviation. To determine the deviation the median of all 
% negative values can used.
%
% To include changes of athmospheric pressure with time, only portions of
% the dive data will be processed at a time. The variable timeframe defines
% the length of every portion in seconds. Ideally this timeframe should
% at least span the longest dive length expected (i.e. 300s for a maximum of
% 5min dives etc.) 
depth=baselinetracking(depth);
N=length(depth);

cdr=strsplit(currpath,'\');
cdr=cdr{length(cdr)};
disp(['Now processing: ' cdr]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%         COMBINE DIVE & GPS DATA         %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for native GPS data in Sensor matrix
gd = find(~isnan(sLat));

if(gd>0)
    disp('Native GPS in sensor file found.');
    GPSData = [timestamp(:), depth(:), sLat(:) sLon(:) , temperature(:)];
else
    % Redundant when using AxyTreks
    try
        disp(['Combining GPS & dive data']);
        % write trackdata matrix into a single variable:
        % DATE/TIME  DEPTH  LAT  LON
        j=1;
        k=1;
        combined=0;
        
        % initializing to expected size boosts performance
        GPSData = [timestamp(:), depth(:), ones(N,2)*nan , temperature(:)];
        
        for j=1:length(timestampGPS)
            row = find(timestamp==timestampGPS(j));
            if ~isempty(row)
                GPSData(row,3:4)=[latGPS(j) lonGPS(j)];
                combined=combined+1;
            end
        end
        
        disp([num2str(combined) ' matches found'])
        
    catch
        disp(['sensor & GPS data not combined (how without valid GPS data???)']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%         SET ANALYSIS PARAMETERS         %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tracer = 1;  % indicates which line in our datafile we are in
eventtimer = 0; % resetable counter for event anlysis
diveno = 0; % number of dive event
eventno = 0; % number of analysed dive event
surfstart = 0; % startindex of surfacetime interval 
events = [];
surftime=[];
divedist = NaN;

% event threshold (seconds): minimum length of event to qualify as dive
eventthres = 5;
if eventthres < interval; eventthres = interval; end

% % depth threshold (m): when can a pressure change considered to be a dive
% if findstr('0390', name)>0;
%     depththres = 0.6; % Wildlife Computers
% else
    depththres = 0.1; % Earth&Ocean
% end

% Dive threshold (m): which dives do we want to analyse in detail (ie.
% wiggles, vertical veolcity etc.)
divethres = 0.5;

% if a bird is at the surface longer than the onlandtreshold (s) we assume
% that it is on land rather then out at sea (longest surface interval
% recorded for Snares Crested Penguins during an overnight trip was 7500s)
% For Yellow-eyed penguins 5.5 hours are a good overnight threshold
onlandthreshold =18000; %5 hours

tripno=1;

wdata = [];
ddata = [];
ndata = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
%         FIND A DIVE EVENT         %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



fprintf(1,'Analysing... ');


while tracer<length(depth)   
    % find start of a dive
    while depth(tracer)<depththres
        tracer=tracer+1;
        if tracer>=length(depth); break; end
    end    
    % find end of a dive
    while depth(eventtimer+tracer)>=depththres  
        eventtimer = eventtimer+1;
        if eventtimer+tracer>length(depth); break; end
    end    
   % do we really have a dive or is it surface oscillation?
    if eventtimer*interval<eventthres
       % no, too short, not  a dive, carry on
       eventtimer=0;
       tracer=tracer+1;
    else   
        %yep, dive, should we analyse it?

        diveno = diveno+1;
        evtstart = tracer-1;
        if(evtstart<1) 
            evtstart=1; 
        end
        evtend = eventtimer+tracer;
        if evtend>length(depth); break; end
        [evtmax, maxindex] = max(depth(evtstart:evtend));
        if surfstart~=0
            surftime=datevec(datenum(timestamp(evtstart))-datenum(timestamp(surfstart)));
            surftime = (surftime(3)*24*60*60)+(surftime(4)*60*60)+(surftime(5)*60)+(surftime(6));
            % bird at surface or on land? (surfacetime longer than
            % on-land-trheshold?)
            if surftime>=onlandthreshold
                diveno = 0;
                eventno = 0;
                tripno = tripno+1;
                surftime = NaN;
            end
            % check if we got a GPS fix during the surface period
            one=find(GPSData(:,1)>timestamp(surfstart));
            two=find(GPSData(:,1)<timestamp(evtstart));
            GPSival=intersect(one,two);
            
            GPSfix=find(~isnan(GPSData(GPSival,3)));
            fix=GPSData(GPSival(GPSfix),3:4);
           
            %
            if ~isempty(fix)           % START if~isempty(fix)
                temp=size(fix);
                if temp(1)<2
                    first_lat=fix(1);
                    first_lon=fix(2);
                    last_lat=NaN;
                    last_lon=NaN;
                    surfdist=0;
                    lastfix=fix;
                    nooffixes=1;
                else
                    if temp(1)==2
                        % This is the 'cheap' distance calculation (i.e. direct
                        % line between first and last fix during surface
                        % interval
                        surfdist=pos2dist(fix(1,1),fix(1,2),fix(temp(1),1),fix(temp(1),2),1)*1000;
                    else
                        % Now this is the 'good' distance because this time we
                        % caluclate the cumulative distances between all 
                        % consecutive fixes. However, this calculation
                        % should become irrelevant once the GPS loggers
                        % firmware works as its supposed to (i.e. only one fix
                        % between dives).
                        cdist=0;
                        for fx = 1:temp(1)-1
                            surfdist=pos2dist(fix(fx,1),fix(fx,2),fix(fx+1,1),fix(fx+1,2),1)*1000;
                            cdist=cdist+surfdist;
                        end
                        surfdist=cdist;
                        clear cdist fx
                    end                   
                    first_lat = fix(1,1);
                    first_lon = fix(1,2);
                    last_lat=fix(temp(1),1);
                    last_lon=fix(temp(1),2);
                    nooffixes=temp(1);
                end
                
            else                       % if~isempty(fix)
                nooffixes=0;
                surfdist=NaN;
                first_lat=NaN;
                first_lon=NaN;
                last_lat=NaN;
                last_lon=NaN;
            end                        % END if~isempty(fix)

        else 
                nooffixes=NaN;
                surfdist=NaN;
                first_lat=NaN;
                first_lon=NaN;
                last_lat=NaN;
                last_lon=NaN;
        end

        if isempty(surftime)
                surftime=NaN;
        end
        if evtmax>=divethres    % Analyse dive!
            eventno=eventno+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                    %
            %         DETERMINE ASCEND AND DESCEND PHASES        %
            %                                                    %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            divedepth = depth(evtstart:evtend);  % define pressure datarange of dive event
            divetimer = timestamp(evtstart:evtend); % define time datarange of dive event
            startbottom = 1; 
            endbottom = length(divedepth);
            descvelo = [];
            ascvelo = [];
            
            % Novel approach: define bottom phase as amount of time between 90% and 100% of max depth)    
%             mx = max(divedepth);
%             bd = find(divedepth>mx*0.9);
%             startbottom = bd(1)-1;
%             endbottom = bd(end)+1;
            
           
            [dd,ddd]=bottomtrap(divedepth);
            startbottom=ddd(1);
    
            endbottom=ddd(2);   
            

     
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                %
            %         WIGGLE ANALYSIS        %
            %                                %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if(endbottom>startbottom) 
            
               wiggledata = divedepth(startbottom:endbottom);
               wigglesummary = [];
               wiggle = 0; 
               noofwiggles = 0; 
               wigglesdetail = [];
               prevdir = 0;  % previous direction of wiggle: -1 = down, +1 = up
               wiggletracer = 1;

               % Wiggleanalysis Part 1: determine amplitudes and velocities of
               % consecutive data points and index them according to previous
               % idrection, i.e. if previous amplitude has the same direction
               % as the current one, we're still in the same wiggle (no
               % direction change) therfore the amp and velo should have the
               % same index as the previous value. 
               % In other words, we end up with a matrix that contains all
               % amplitude and velocitiy values of the entire bottom phase but
               % with a wigglesindex in the first column

               while wiggletracer+1<=length(wiggledata)    
                   amp = wiggledata(wiggletracer+1)-wiggledata(wiggletracer);
                   velo = amp/interval;

                   if wiggledata(wiggletracer)>wiggledata(wiggletracer+1)
                       if prevdir>0
                           wigglesdetail = [ wigglesdetail; wiggle wiggledata(wiggletracer) amp velo ];
                       else
                           wiggle = wiggle+1;
                           wigglesdetail = [ wigglesdetail; wiggle wiggledata(wiggletracer) amp velo];                         
                       end                   
                       prevdir = 1;
                   elseif wiggledata(wiggletracer)<wiggledata(wiggletracer+1)
                       if prevdir<0
                           wigglesdetail = [ wigglesdetail; wiggle wiggledata(wiggletracer) amp velo];
                       else
                           wiggle = wiggle+1;
                           wigglesdetail = [ wigglesdetail; wiggle wiggledata(wiggletracer) amp velo];                         
                       end                   
                       prevdir = -1;
                   else  
                   end
                   wiggletracer=wiggletracer+1;     
               end

               % Wiggleanalysis Part 2: now we have to draw the wiggle
               % parameters from the wiggle matrix, i.e. how many wiggles
               % occured, what is the duration of every wiggle, what is 
               % the vertical velocity for every wiggle, what is the
               % amplitude of every wiggle.
               % We will write only a summary of the wiggle analysis into the
               % main result file. Detailed wiggle summary will be written into
               % a seperate wiggle-file.

               if isempty(wigglesdetail)<1
                 noofwiggles = max(wigglesdetail(:,1)); 
                 % the variable WIGGLESUMMARY will contain the the following values:
                 % INDEX   ONSET DEPTH   LENGTH OF WIGGLE   TOTAL AMPLITUDE OF WIGGLE   MEAN VERTICAL VELOCITY

                 for i = 1:noofwiggles
                       currentwiggle = find(wigglesdetail(:,1)==i);
                       if wigglesdetail(currentwiggle,3)<0
                           onsetdepth = max(wigglesdetail(currentwiggle,2));
                       else
                           onsetdepth = min(wigglesdetail(currentwiggle,2));
                       end

                      wigglesummary = [ wigglesummary; i (length(currentwiggle))*interval onsetdepth sum(wigglesdetail(currentwiggle,3)) mean(wigglesdetail(currentwiggle,4)) std(wigglesdetail(currentwiggle,4))];

                      wiggleline =  [tripno, eventno, timestamp(evtstart), timestamp(evtstart), timestamp(evtend), evtmax, i, wigglesummary(i,2), ...
                            wigglesummary(i,3), wigglesummary(i,4), wigglesummary(i,5), wigglesummary(i,6)];


                       wdata = [wdata; wiggleline ];
    %                    fprintf(fidwiggles,newline);
    %                    fprintf(fidwiggles,'\n');
                 end
               end
            end
           
           
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            % Write analysis details to files %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            TripNumber=  [TripNumber; tripno];
            EventNo=     [EventNo; eventno];
            Date=        [Date; datestr(timestamp(evtstart),20)];
            SurfaceTime= [SurfaceTime; surftime];
            NoFixes=     [NoFixes; nooffixes];
            FirstLat=    [FirstLat; first_lat];
            FirstLon=    [FirstLon; first_lon];
            LastLat=     [LastLat; last_lat];
            LastLon=     [LastLon; last_lon];
            SurfaceDist= [SurfaceDist; surfdist];
            SurfaceVelo= [SurfaceVelo; surfdist/surftime];
            EvtOnset=    [EvtOnset; datestr(timestamp(evtstart),13)];
            EvtEnd=      [EvtEnd; datestr(timestamp(evtend),13)];
            DiveTime=    [DiveTime; eventtimer*interval];
            DiveDistance=[DiveDistance; NaN];
            DiveSpeed=   [DiveSpeed; NaN];
            EvtMaxDepth= [EvtMaxDepth; evtmax];
            DescDur=     [DescDur; startbottom*interval];
            if(startbottom>0)
                DescVelo=    [DescVelo; divedepth(startbottom)/(startbottom*interval)];
                StartBot=    [StartBot; divedepth(startbottom)];
            else
                DescVelo=    [DescVelo; 0];
                StartBot=    [StartBot; 0];
            end
            
            NoWiggles=   [NoWiggles; noofwiggles];
            BottomTime=  [BottomTime; (endbottom-startbottom)*interval];
            if(endbottom>0)
                BottomEnd=   [BottomEnd; divedepth(endbottom)];
            else
                BottomEnd=   [BottomEnd; 0];
            end
            AscDur=      [AscDur; (length(divedepth)-endbottom)*interval];
            if(endbottom>0)
                AscVelo=     [AscVelo; divedepth(endbottom)/((length(divedepth)-endbottom)*interval)];
            else
                AscVelo=     [AscVelo; 0];
            end

            
                  if ~isempty(wigglesummary)
                MeanWigAmp=  [MeanWigAmp; mean(wigglesummary(:,4))];
                MinWigAmp=   [MinWigAmp; min(wigglesummary(:,4))];
                MaxWigAmp=   [MaxWigAmp; max(wigglesummary(:,4))];
                MeanWigVelo= [MeanWigVelo; mean(wigglesummary(:,5))];
                MinWigVelo=  [MinWigVelo; min(wigglesummary(:,5))];
                MaxWigVelo=  [MaxWigVelo; max(wigglesummary(:,5))];
            else
                MeanWigAmp=  [MeanWigAmp; NaN];
                MinWigAmp=   [MinWigAmp; NaN];
                MaxWigAmp=   [MaxWigAmp; NaN];
                MeanWigVelo= [MeanWigVelo; NaN];
                MinWigVelo=  [MinWigVelo; NaN];
                MaxWigVelo=  [MaxWigVelo; NaN];
            end
            
          

            events = [ events; eventno evtstart evtstart+startbottom-1 evtstart+endbottom-1 evtend ] ;
            surfstart = evtend;
            tracer=evtend+1;
            eventtimer=0;
        else                % Do not analyse dive, just give out basic parameters

            noiseline =  [diveno surftime timestamp(evtstart) timestamp(evtend), eventtimer*interval, evtmax];
            ndata = [ndata; noiseline];

            eventtimer=0;
            tracer=evtend+1;
        end
    end
    done = (tracer-1)/length(depth);

end

ddata= struct('TripNumber',TripNumber,'EventNo',EventNo,'Date',Date,...
    'SurfaceTime',SurfaceTime,'NoFixes',NoFixes,'FirstLat',FirstLat,...
    'FirstLon',FirstLon,'LastLat',LastLat,'LastLon',LastLon,...
    'SurfaceDist',SurfaceDist,'SurfaceVelo',SurfaceVelo,'EvtOnset',EvtOnset,...
    'EvtEnd',EvtEnd,'DiveTime',DiveTime,'DiveDistance',DiveDistance,...
    'DiveSpeed',DiveSpeed,'EvtMaxDepth',EvtMaxDepth,'DescDur',DescDur,...
    'DescVelo',DescVelo,'StartBot',StartBot,'NoWiggles',NoWiggles,...
    'MeanWigAmp',MeanWigAmp,'MinWigAmp',MinWigAmp,'MaxWigAmp',MaxWigAmp,...
    'MeanWigVelo',MeanWigVelo,'MinWigVelo',MinWigVelo,'MaxWigVelo',MaxWigVelo,...
    'BottomTime',BottomTime,'BottomEnd',BottomEnd,'AscDur',AscDur,...
    'AscVelo',AscVelo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                       %
% CALCULATE DIVING & FORAGING EFFICIENY %
%                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tn=ddata.TripNumber;
st=ddata.SurfaceTime;
dt=ddata.DiveTime;
bt=ddata.BottomTime;

efficiency(1:length(tn),1:2)=NaN;


for i=1:length(tn)
  if ~isnan(st(i))
    efficiency(i,1)=dt(i)/(st(i)+dt(i));
    efficiency(i,2)=bt(i)/(st(i)+dt(i));
  end
end

ddata.DivingEfficiency = efficiency(:,1);
ddata.ForagingEfficiency = efficiency(:,2);

clear tn st dt bt efficiency


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                            %
% CALCULATE DIVE DISTANCE FROM GPS POSITIONS %
%                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Re-caluclate DiveDistance
    temp=[ddata.NoFixes ddata.FirstLat...
        ddata.FirstLon,ddata.LastLat...
       ddata.LastLon];
    DiveDistance(1:length(temp),1)=NaN;
    for i =1:size(temp,1)-1
        if temp(i,1)==1
            DiveDistance(i,1)=pos2dist(temp(i,2),temp(i,3),temp(i+1,2),temp(i+1,3),1)*1000;
        elseif temp(i,1)>1
            DiveDistance(i,1)=pos2dist(temp(i,4),temp(i,5),temp(i+1,2),temp(i+1,3),1)*1000;
        end
    end
    DiveDistance=[DiveDistance; NaN];
    ddata.DiveDistance=DiveDistance;
    clear temp i DiveDistance;
    
    
%%%%%%%%%%%%%%%%%%%%
%                  %
% FIX TRIP NUMBERS %
%                  %
%%%%%%%%%%%%%%%%%%%%

temp = unique(ddata.TripNumber);
ntn = 1:1:length(temp);

ntn=ntn';

for i=1:length(temp)
    j=find(ddata.TripNumber==temp(i));
    ddata.TripNumber(j)=ntn(i);
end

% Throw out bogus trips


clear temp ntn i j;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     %
% Calculate homerange, trip distance & segment speeds %
%                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prioritize LastLat/lons over FirstLat/lons
trips = unique(ddata.TripNumber);
fixes = [ddata.FirstLat ddata.FirstLon];
fixesB = [ddata.LastLat ddata.LastLon];
bf = find(~isnan(fixesB(:,1)));
fixes(bf,:)=fixesB(bf,:);

clear fixesB bf;




% a) Calcuate distance from home
homerange(1:length(fixes),1) = NaN;
for i=1:length(fixes)
    if ~isnan(fixes(i,1))
        homerange(i,1)=pos2dist(home(1),home(2),fixes(i,1),fixes(i,2),1);
        
    else
        homerange(i,1)=NaN;
    end
end
ddata.HomeRange=homerange;
clear homerange;

% Important: differentiate between trips
cumuldist(1:length(fixes),1) = NaN;
for t=1:length(trips)
    % b) calculate cumulative distance
    tripInd = find(ddata.TripNumber==trips(t));
    tripfixes = fixes(tripInd,:); 
    if(size(tripfixes,1)<5) 
        continue %skip super short trips
    end
    prevLat=home(1);
    prevLon=home(2);
    currdist = 0;
    for i=1:length(tripfixes)
        if ~isnan(tripfixes(i,1))
            cumuldist(tripInd(i),1) = currdist+pos2dist(prevLat,prevLon,tripfixes(i,1),tripfixes(i,2),1);
            currdist=cumuldist(tripInd(i),1);
            prevLat=tripfixes(i,1);
            prevLon=tripfixes(i,2);
        end
    end
    cumuldist(tripInd(end),1) = currdist+pos2dist(prevLat,prevLon,home(1),home(2),1);
end

ddata.CumulativeDistance=cumuldist;
clear currdist cumuldist;

% c) calculate segment speed

segmentspeed(1:length(fixes),1) = NaN;

temp_dmy = datevec(ddata.Date,'dd/mm/yy');
temp_hms = datevec(ddata.EvtOnset);
tempTS = datenum([temp_dmy(:,1:3) temp_hms(:,4:6)]);

fixes=[fixes tempTS];
clear temp_dmy temp_hms tempTS;

prevLat=home(1);
prevLon=home(2);
prevTime=fixes(1,3);

for i=1:length(fixes)
    if ~isnan(fixes(i,1))
        dist = pos2dist(prevLat,prevLon,fixes(i,1),fixes(i,2),1)*1000; %dist in metres
        timediff = (fixes(i,3)-prevTime)*86400; %time in seconds
        segmentspeed(i,1)=dist/timediff;

        prevLat=fixes(i,1);
        prevLon=fixes(i,2);
        prevTime=fixes(i,3);
    end
end
dist = pos2dist(prevLat,prevLon,home(1),home(2),1)*1000; %dist in metres
timediff = (fixes(i,3)-prevTime)*86400; %time in seconds
segmentspeed(i,1)=dist/timediff;
ddata.SegmentSpeed=segmentspeed;
clear dist timediff prevTime prevLat prevLon;

    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
% EXTRAPOLATE DIVE POSITIONS %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% variable 'fixes' carries over from above

% paste home location as first and last fix (if needed)
if isnan(fixes(1,1))
    fixes(1,1)=home(1);
    fixes(1,2)=home(2);
end

if isnan(fixes(end,1))
    fixes(end,1)=home(1);
    fixes(end,2)=home(2);
end

fixIdx = find(~isnan(fixes(:,1)));

expolFixes(1:length(fixes),1:4) = NaN;

for i=1:length(fixIdx)-1
    a=fixIdx(i);
    b=fixIdx(i+1);
    
    segLatDiff=fixes(b,1)-fixes(a,1);
    segLonDiff=fixes(b,2)-fixes(a,2);
    segTimDiff=fixes(b,3)-fixes(a,3);
    
    %temporal gap size
    gap=segTimDiff*86400;
    
    j=a+1;
    while j<b
    %for j=(a+1):(b-1)
       fixTime = fixes(j,3);
       relTime = (fixTime-fixes(a,3))/segTimDiff;
       relLatInc = segLatDiff*relTime;
       relLonInc = segLonDiff*relTime;
       
       expolFixes(j,1)=fixes(a,1)+relLatInc;
       expolFixes(j,2)=fixes(a,2)+relLonInc;
       
       %calculate TTF
       tempTTFa = fixTime-fixes(a,3);
       tempTTFb = fixes(b,3)-fixTime;
       if tempTTFb>tempTTFa
           expolFixes(j,3)=tempTTFa*86400;
       else
           expolFixes(j,3)=tempTTFb*86400;
       end
       
       expolFixes(j,4)=gap;
       
       j=j+1;
    end 
    
end


ddata.Lat = fixes(:,1);
ddata.Lon = fixes(:,2);
ddata.ExpolLat = expolFixes(:,1);
ddata.ExpolLon = expolFixes(:,2);
ddata.ExpolTTF = expolFixes(:,3);
ddata.ExpolGap = expolFixes(:,4);





fprintf(1,'...done!\n\n');

fname=name;
