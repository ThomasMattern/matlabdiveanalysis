%THIS SCRIPT ASSUMES THAT THE DATA FILES FOR THE DIFFERENT BIRDS ARE STORED
%IN SEPERATE FOLDER WHICH ARE NAMED ACCORDING TO BIRD ID, DEPLOYMENT DATE
%ETC.
note=0;
% get data dir
path = uigetdir('D:\MATLAB\work\');


home = [];

[indx,~] = listdlg('PromptString',{'Select site:'},...
    'SelectionMode','single','ListString',{'Harrison Cove','Moraine',...
    'East Shelter Is','Seymour Is','Rolla Is','Pigeon House','Groper Is',...
    'Whenua Hou','Jackson Head','Charleston','Te Kahui','Choros A1',...
    'Choros A2','Choros A4','Choros A5','Choros A9','Choros W',...
    'Anchorage Bay East','South Bay West','Stack Bay','Te Kāhui',...
    'El Pedral','Many Is #3','Petrel Is'},'ListSize',[150,120]);
switch indx
    case 1
        home=[-44.622771,167.912524]; %Harrison Cove
    case 2
        home=[-44.604704,167.810033]; %Moraine
    case 3
        home=[-45.271234,166.894132]; %East Shelter
    case 4
        home=[-45.307588,167.006702]; %Seymour Is
    case 5
        home=[-45.441135,167.132435]; %Rolla Is
    case 6
        home=[-47.226825,167.653951]; %Pigeon House
    case 7
        home=[-46.953055,168.144873]; %Groper Is
    case 8  
        home=[-46.765129,167.647371]; %Whenua Hou
    case 9  
        home=[-43.961788,168.613784]; %Jackson Head
    case 10  
        home=[-41.893092,171.446476]; %Charleston
    case 11  
        home=[-39.05924089608756, 174.04805297603204]; %Te Kahui
    case 12  
        home=[-29.259357, -71.535535]; %Choros A1
    case 13  
        home=[-29.257004, -71.536251]; %Choros A2   
    case 14  
        home=[-29.255450, -71.544577]; %Choros A4   
    case 15  
        home=[-29.265617, -71.538186]; %Choros A5 
    case 16  
        home=[-29.250559, -71.536961]; %Choros A9 
    case 17  
        home=[-29.263734, -71.544799]; %Choros W   
    case 18  
        home=[-49.668056, 178.807491]; %Anchorage Bay East 
    case 19  
        home=[-49.707234, 178.748319]; %South Bay West   
    case 20  
        home=[-49.694825, 178.740793]; %Stack Bay       
    case 21  
        home=[-39.059280, 174.048080]; %Te Kahui   
    case 22  
        home=[-42.946349, -64.363470]; %El Pedral
    case 23
        home=[-45.773936,166.512926]; %Many Is #3
    case 24
        home=[-45.741263,166.519396]; % Petrel Is
end
clear indx tf;

%home=[-45.863520,170.741134];
        
% dlg = 'Enter comma-separated home coordinates (i.e. lat,lon):';
% while(length(home)~=2)
%     homeINP = inputdlg(dlg,...
%                      'Home coordinates', [1 60]);
%     home = str2num(homeINP{:});
%     if(length(home)~=2)
%         dlg = ['INPUT ERROR. Format must be lat,lon'];
%     end
% end
% clear dlg homeINP;

% get contents of dir
dircont = dir(path);
% delete . and .. from the beginning of the directory listing
dircont(1:2)=[];
% delete non-directory entries from dir listing
rgi=1;
kill=[];
while rgi <= length(dircont)
    if dircont(rgi).isdir==0
        kill=[kill; rgi];
    end
    rgi=rgi+1;
end
dircont(kill)=[];
clear kill;


DiveData(:).birdID=[];
DiveData(:).dive=[];
DiveData(:).sensor=[];
DiveData(:).gps=[];


% extract the SENSOR/GPS files consecutively from the directories
rgi = 1;
while rgi <= length(dircont)
    name=dircont(rgi).name;
    
    currpath=[path '\' dircont(rgi).name];
    currdir=dir(currpath);
    currdir(1:2)=[];
    
    rgj = 1;
    while rgj <= length(currdir)
    
        fname = currdir(rgj).name;
        nl = length(currdir(rgj).name);
        %type = currdir(rgj).name(nl-6:nl-4);
        type = fname(end-6:end-4);

        if(type=='ive')
            try
               [sday,smonth,syear,shour,sminute,ssecond,temperature,depth,sLat,sLon]...
                   = textread([currpath '\' fname],...
                   '%d-%d-%d,%d:%d:%d,%f,%f,%f,%f%*[^\n]',...
                   'headerlines',1);
            catch
                [sday,smonth,syear,shour,sminute,ssecond,temperature,depth]...
               = textread([currpath '\' fname],...
               '%d-%d-%d,%d:%d:%d,%f,%f%*[^\n]',...
               'headerlines',1);
                clear sLat sLon;
                sLat(1:length(sday),1) = NaN;
                sLon(1:length(sday),1) = NaN;
            end
           % check depth values. If in bar calculate average minimum and
           % use to convert to depth [m]
           if(max(depth)>500)
               %mins = find(depth<1012);
               %surfpres = mean(depth(mins));
               %depth = depth/surfpres;
               depth = (depth-1012)/100;
           end    
           
           SENSOR = [sday smonth syear shour sminute ssecond temperature depth sLat sLon];
           clear sday smonth syear shour sminute ssecond temperature depth sLat sLon;
        elseif (type=='gps')
            if isfile([currpath '\' fname])
                try
                   [gday,gmonth,gyear,ghour,gminute,gsecond,lat,lon]...
                       = textread([currpath '\' fname],...
                       '%2d.%2d.%4d,%2d:%2d:%2d\t%f\t%f%*[^\n]');
                catch
                    try
                        [gday,gmonth,gyear,ghour,gminute,gsecond,lat,lon]...
                              = textread([currpath '\' fname],...
                            '%2d/%2d/%4d,%2d:%2d:%2d\t%f\t%f%*[^\n]');  
                    catch
                        [gday,gmonth,gyear,ghour,gminute,gsecond,lat,lon]...
                             = textread([currpath '\' fname],...
                             '%2d/%2d/%4d\t%2d:%2d:%2d\t%f\t%f%*[^\n]'); 
                    end
                end
                GPS = [gday gmonth gyear ghour gminute gsecond lat lon];
                 clear gday gmonth gyear ghour gminute gsecond lat lon;
            end
        end
        rgj = rgj+1;
    end
    


    
    % HERE GOES THE ANALYSIS STUFF

    
    dive_analysis
    
    DiveData(rgi).birdID=[name];
    DiveData(rgi).dive=ddata;
    DiveData(rgi).sensor=SENSOR;
    if exist('GPS') 
        DiveData(rgi).gps=GPS;
    end
    %struct('birdID',name,'dive',ddata,'wiggles',wdata,'noise',ndata,'sensor',SENSOR,'gps',GPS);
    rgi=rgi+1;
   
    %clear ndata ddata name wdata
end
%dive_stats_YEP

write_data

% OPTIONAL: write DiveData into MS Office 2003 Excel XML speradsheet
% (without the wiggle data, because that is simply too much for Bill's
% sorftware)
%resname=[fname '_DiveAnalysis_'];
%wrt_YEP
%write_excel_xml(resname,DiveData)