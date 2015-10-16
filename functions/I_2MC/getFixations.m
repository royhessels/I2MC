function [cutoff,starttime,endtime,fixdur,xmedian,ymedian,flankdataloss,fracinterped,RMSxy,BCEA,fixRangeX,fixRangeY] = getFixations(finalweights,missing,xpos,ypos,timestamp,cutoffstd,pixperdeg,maxMergeDist,maxMergeTime,minFixDur)
% determine fixations based on finalweights from 2-means clustering

% Roy Hessels - 2014

%% input:

% finalweights              = 2-means clustering weighting
% missing                   = boolean vector indicating which samples are
%                             missing (originally, before interpolation!)
% xpos, ypos                = horizontal & vertical coordinates from ET
% timestamp                 = timestamp from ET (should be in ms!)
% cutoffstd                 = number of std above mean clustering-weight to
%                               use as fixation cutoff
% pixperdeg                 = used for transforming pixels on screen to
%                               degrees
% maxMergeDist              = maximum Euclidean distance in pixels between fixations for merging
% maxMergeTime              = maximum time in ms between fixations for merging

%% output:

% cutoff                    = cutoff used for fixation detection
% starttime                 = vector with fixation start times
% endtime                   = vector with fixation end times
% fixdur                    = vector with fixation durations
% xmedian                   = vector with fixation median horizontal
%                               position (one value for each fixation in
%                               trial)
% ymedian                   = vector with fixation median vertical
%                               position (one value for each fixation in
%                               trial)
% flankdataloss             = boolean with 1 for when fixation is flanked
%                               by data loss, 0 if not flanked by data loss
% RMSxy                     = RMS of fixation (precision)
% BCEA                      = BCEA of fixation (precision)
% fixRangeX                 = max(xpos) - min(xpos) of fixation
% fixRangeY                 = max(ypos) - min(ypos) of fixation

%% output of previous version (before 11/03/2015):

% fixstart                  = vector with fixation start indices
% fixend                    = vector with fixation end indices

%% first determine cutoff for finalweights
cutoff = nanmean(finalweights) + cutoffstd*nanstd(finalweights);

% get boolean of fixations
fixbool = finalweights < cutoff;

% get indices of where fixations start and end
[fixstart,fixend] = bool2bounds(fixbool);

% get start time, end time, and fix duration
starttime   = timestamp(fixstart);
endtime     = timestamp(fixend);

%% loop over all fixations in trial, see if should be merged
for p=length(starttime):-1:2
    % get median coordinates of fixation
    xmedThis = median(xpos(fixstart(p):fixend(p)));
    ymedThis = median(ypos(fixstart(p):fixend(p)));
    xmedPrev = median(xpos(fixstart(p-1):fixend(p-1)));
    ymedPrev = median(ypos(fixstart(p-1):fixend(p-1)));
    % check if fixations close enough in time and space and thus qualify
    % for merging
    if starttime(p)-endtime(p-1)<maxMergeTime && hypot(xmedThis-xmedPrev,ymedThis-ymedPrev) < maxMergeDist
        % merge
        fixend(p-1) = fixend(p);
        endtime(p-1)= endtime(p);
        % delete merged fixation
        fixstart(p) = [];
        fixend(p)   = [];
        starttime(p)= [];
        endtime(p)  = [];
    end
end

%%  beginning and end of fixation must be real data, not interpolated. If 
% interpolated, those bit(s) at the edge(s) is excluded from the fixation
% first throw out fixations that are all missing/interpolated
for p=length(starttime):-1:1
    if all(missing(fixstart(p):fixend(p)))
        fixstart(p) = [];
        fixend(p)   = [];
        starttime(p)= [];
        endtime(p)  = [];
    end
end
for p=1:length(starttime)
    if missing(fixstart(p))
        fixstart(p) = fixstart(p)+find(~missing(fixstart(p):   fixend(p)  ),1)-1;
        starttime(p)= timestamp(fixstart(p));
    end
    if missing(fixend(p))
        fixend(p)   = fixend(p)  -find(~missing(fixend(p)  :-1:fixstart(p)),1)+1;
        endtime(p)  = timestamp(fixend(p));
    end
end


%% check if any too short
fixdur      = endtime-starttime;
qTooShort   = fixdur<minFixDur;
fixstart(qTooShort) = [];
fixend(qTooShort)   = [];
starttime(qTooShort)= [];
endtime(qTooShort)  = [];
fixdur(qTooShort)   = [];

%% process fixations, get other info about them
[   xmedian,...         % vectors for medians
    ymedian,...
    flankdataloss,...   % vectors for whether fixation is flanked by data loss
    fracinterped,...
    RMSxy,...           % vectors for precision measures
    BCEA,...
    fixRangeX,...
    fixRangeY]      = deal(zeros(size(fixstart)));

for a = 1:length(fixstart)
    % get data during fixation
    xposf = xpos(fixstart(a):fixend(a));
    yposf = ypos(fixstart(a):fixend(a));
    % for all calculations below we'll only use data that is not
    % interpolated, so only real data
    qMiss = missing(fixstart(a):fixend(a));
    
    % get median coordinates of fixation
    xmedian(a) = median(xposf(~qMiss));
    ymedian(a) = median(yposf(~qMiss));
    
    % determine whether fixation is flanked by period of data loss
    flankdataloss(a) = (fixstart(a)>1 && missing(fixstart(a)-1)) || (fixend(a)<length(xpos) && missing(fixend(a)+1));
    
    % fraction of data loss/interpolated data during fixation
    % assumption: missing data during fixation has been interpolated and is
    % thus no longer missing
    % assert to check that assumption:
    assert(all(~isnan(xposf(qMiss))),'assumption failed')
    fracinterped(a)  = sum(qMiss)./(fixend(a)-fixstart(a)+1);
    
    % calculate RMS
    % since its done with diff, don't just exclude missing and treat
    % resulting as one continuous vector. replace missing with nan first,
    % use left-over values
    xdif = subsasgn(xposf,substruct('()',{qMiss}),nan);
    ydif = subsasgn(yposf,substruct('()',{qMiss}),nan);
    xdif = diff(xdif).^2; xdif(isnan(xdif)) = [];
    ydif = diff(ydif).^2; ydif(isnan(ydif)) = [];
    c    = xdif + ydif; % 2D sample-to-sample displacement value in pixels
    visuelehoeksq = c./pixperdeg^2; % value in degrees visual angle
    RMSxy(a) = sqrt(mean(visuelehoeksq));
    
    % calculate BCEA (Crossland and Rubin 2002 Optometry and Vision Science)
    stdx = std(xposf(~qMiss))/pixperdeg; % value in degrees visual angle
    stdy = std(yposf(~qMiss))/pixperdeg; % value in degrees visual angle
    xx   = corrcoef(xposf(~qMiss),yposf(~qMiss));
    rho  = xx(1,2);
    P    = 0.68; % cumulative probability of area under the multivariate normal
    k    = log(1/(1-P));
    
    BCEA(a) = 2*k*pi*stdx*stdy*sqrt(1-rho.^2);
    
    % calculate max-min of fixation
    fixRangeX(a) = (max(xposf(~qMiss)) - min(xposf(~qMiss)))/pixperdeg;
    fixRangeY(a) = (max(yposf(~qMiss)) - min(yposf(~qMiss)))/pixperdeg;
end

return