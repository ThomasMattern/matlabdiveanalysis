function resdepth = baselinetracking(depth)
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

errorrange = 0.5; 
timeframe=600; % track baseline in this timeframe throughout the data

if length(depth)<timeframe
    timeframe=length(depth);
end
bits = length(depth)/timeframe;
if bits<round(bits)
    bits = round(bits)-1;
else
    bits=round(bits);
end
i = 1;
currentbit=1;
depthmode = [];
N=length(depth);

while i < length(depth)
    block = depth(i:i+(timeframe-1));
    blockvalues = find(block<=errorrange);
    depthmode = median(block(blockvalues));
    if isempty(depthmode)
        depthmode=median(find(block<=min(block)+5));
    end

    if depthmode<0
        block=block+depthmode;
        block(find(block<=0))=0;
        depth(i:i+timeframe-1)=block;
          else
        block=block-depthmode;
        block(find(block<=0))=0;
        depth(i:i+timeframe-1)=block;
     end    

    i=i+timeframe; 
    currentbit = currentbit+1;
    if currentbit>bits 
        timeframe=length(depth)-i; 
    end  
end

killer=find(depth<0);
if ~isempty(killer)
    depth(killer)=0;
end

resdepth = depth;
