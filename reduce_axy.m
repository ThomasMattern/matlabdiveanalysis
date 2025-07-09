% Script to prepare data recorded with AxyTrek GPS dive loggers for dive analyis.
% The script reduces the dive/accelerometer data file so that only the rows with
% dive data remain. 


[filename, pathname] = uigetfile('.csv', 'Select raw AxyDepth dive file: ');
rawdive = [pathname filename];
clear pathname filename;

datetime.setDefaultFormats('default','dd/MM/yyyy');

ds = tabularTextDatastore(rawdive,'TextscanFormats',{'%q',...
    '%{dd/MM/uuuu}D','%T','%f','%f','%f','%q','%f','%f','%f',...
    '%f','%f','%f','%f','%f','%f'});
reset(ds)
X = [];
dlm = [];
 disp('Simmering data down. Please be patient.');
while hasdata(ds)
      T = read(ds);
      rows = find(~isnan(T.Depth));
      X = [X; T(rows,:)];
end

if isempty(dlm)
    datetxt=char(X.Date(1));
    dlm=regexp(datetxt,'[^a-z0-9_]');
    dlm=datetxt(dlm(1));
end

bogus = find(X.Date=="01/01/0001");
 X(bogus,:)=[];

% Determine date arrangement
datecompsA = split(char(X.Date(1)),dlm);
datecompsB = split(char(X.Date(end)),dlm);

if(str2num(datecompsA{1})>12 || str2num(datecompsB{1})>12) 
    % first number is not the month
    dtf=['dd' dlm 'mm' dlm 'yyyy'];
elseif(str2num(datecompsA{2})>12 || str2num(datecompsB{2})>12)     
    % second number is not the month
    dtf=['mm' dlm 'dd' dlm 'yyyy'];
elseif(str2num(datecompsB{1})-str2num(datecompsA{1})==0)
    % this only works if all data recorded in same month
    dtf=['mm' dlm 'mm' dlm 'yyyy']; 
elseif(str2num(datecompsB{2})-str2num(datecompsA{2})==0)
    dtf=['dd' dlm 'mm' dlm 'yyyy']; 
else
    dtf = questdlg(['Date in the first line is ' X.Date(1) '. Which of the following date formats is correct?'], ...
	'Clarify date format', ...
	['dd' dlm 'mm' dlm 'yyyy'],['mm' dlm 'dd' dlm 'yyyy'],['dd' dlm 'mm' dlm 'yyyy']);
end


%dtf=['dd' dlm 'mm' dlm 'yyyy'];
dt=datevec(char(X.Date),dtf);
tm=datevec(X.Time,'HH:MM:ss');
dttm = [dt(:,1:3) tm(:,4:6)];
clear dt tm
xaxis = datetime(dttm,'timezone','UTC');
xaxis.TimeZone="America/Santiago";
%xaxis.TimeZone="Pacific/Auckland";
yaxis = -X.Depth;
taxis =X.Temp___C_;
%plot(xaxis,yaxis)
%xtickformat('HH:mm');

ax1 = subplot('Position', [0.1, 0.75, 0.8, 0.18]); % [left, bottom, width, height]
plot(xaxis,taxis)
set(gca, 'XTickLabel', []);
ax2 = subplot('Position', [0.1, 0.1, 0.8, 0.63]); % [left, bottom, width, height]
plot(xaxis,yaxis)
xtickformat('HH:mm');
linkaxes([ax1, ax2], 'x'); % Link only the x-axis

disp('Writing simmered down data to file. Please wait.');
[fiddive,msg]=fopen([rawdive(1:length(rawdive)-4) '_reduced.txt'],'wt');
fprintf(fiddive,[ 'Date', ',',...
                  'Time',',',...
                  'Temperature', ',',...
                  'Depth', ',',...
                  'lat', ',',...
                  'lon','\n']);
        
              
for r=1:size(X,1)
    fprintf(fiddive,[ datestr(dttm(r,:),'dd-mm-yyyy,HH:MM:ss'), ',',...
                      num2str(X.Temp___C_(r)), ',',...
                      num2str(X.Depth(r)), ',',...
                      num2str(X.location_lat(r)), ',',...
                      num2str(X.location_lon(r)), '\n']);
end
fclose all;
clear all;
