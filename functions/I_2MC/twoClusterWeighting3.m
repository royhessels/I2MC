function [finalweights,stopped] = twoClusterWeighting3(xpos,ypos,windowtime,steptime,freq,missingxvalue,missingyvalue,maxerrors)
% Calculates 2-means cluster weighting for eye-tracking data

% Basis final weighting on mean of:
% - weights from data in original sampling frequency
% - weights from data downsampled to 1/3 of original sampling frequency
% - weights from data downsampled to 1/5 of original sampling frequency
% - weights from data downsampled to 1/10 of original sampling frequency
% Example: if SF=300Hz, weights are based on data in 300Hz, 100Hz, 60Hz & 50Hz

% Input:
% xpos,ypos                     = horizontal and vertical coordinates from eye-tracker over which to calculate 2-means clustering
% windowtime                    = time window (s) over which to calculate 2-means clustering (choose value so that max. 1 saccade can occur)
% steptime                      = time window (s) in each iteration. Use zero for sample by sample processing
% freq                          = sampling frequency of data
% missingxvalue                 = missing value for horizontal position
% missingxvalue                 = missing value for vertical position
% maxerrors                     = maximum number of errors allowed in k-means clustering procedure before proceeding to next file

% Output:
% finalweights                  = vector of 2-means clustering weights (one weight for each sample), the higher, the more likely a saccade happened
% stopped                       = whether maxerrors was reached or not

% Roy Hessels - 2014

% calculate number of samples of the moving window
nrsamples =       round(windowtime/(1/freq));
stepsize  = max(1,round(  steptime/(1/freq)));

% create empty weights vectors
totalweights = zeros([length(xpos) 1]);
nrtests = zeros([length(xpos) 1]);

% stopped is always zero, unless maxiterations is exceeded. this
% indicates that file could not be analysed after trying for x iterations
stopped = 0;
counterrors = 0;

% filter signal. Decimate first does a chebychev filter as specified below
rip = .05;	% passband ripple in dB
[b3,a3] = cheby1(8, rip, .8/3);
[b5,a5] = cheby1(8, rip, .8/5);
[b10,a10] = cheby1(8, rip, .8/10);
idx_3 = fliplr(nrsamples-1:-3:1);
idx_5 = fliplr(nrsamples-2:-5:1);
idx_10 = fliplr(nrsamples-5:-10:1);

% see where are missing in this data, for better running over the data
% below.
[on,off] = bool2bounds(xpos==missingxvalue | ypos==missingyvalue);
if ~isempty(on)
    %  merge intervals smaller than nrsamples long
    merge=find(on(2:end)-off(1:end-1)-1<nrsamples);
    for p=fliplr(merge)
        off(p)   = off(p+1);
        off(p+1) = [];
        on (p+1) = [];
    end
    % check if intervals at data start and end are large enough
    if on(1)<nrsamples+1
        % not enough data point before first missing, so exclude them all
        on(1)=1;
    end
    if off(end)>length(xpos)-nrsamples
        % not enough data points after last missing, so exclude them all
        off(end)=length(xpos);
    end
    % start at first non-missing sample if trial starts with missing (or
    % excluded because too short) data
    if on(1)==1
        i=off(1)+1;  % start at first non-missing
    else
        i=1;
    end
else
    i=1;
end
eind = i+nrsamples-1;
while eind<=length(xpos)
    % check if max errors is crossed
    if counterrors > maxerrors
        fprintf('Too many empty clusters encountered, aborting file. \n');
        stopped = 1;
        finalweights = NaN;
        return
    end
    
    % select data portion of nrsamples
    idx = i:eind;
    ll_1 = [xpos(idx) ypos(idx)];
    
    % the filtering done in decimate is done once above. now we simply need
    % to select each nth sample where n is the integer factor by which
    % number of samples is reduced. select samples such that they are till
    % end of window
    ll_3 = filtfilt(b3,a3,ll_1);
    ll_3 = ll_3(idx_3,:);
    ll_5 = filtfilt(b5,a5,ll_1);
    ll_5 = ll_5(idx_5,:);
    ll_10 = filtfilt(b10,a10,ll_1);
    ll_10 = ll_10(idx_10,:);
    
    if 0
        % plot original and down-sampled data
        % note filtering introduced a means shift, but that does not matter
        % as we are interested in switches
        figure(1), clf
        subplot(1,3,1), cla
        hold on
        plot(ll_1(:,1),'b')
        plot(ll_3(SmartVec(1:1:20,3,'flat'),1),'r')
        plot(ll_5(SmartVec(1:1:12,5,'flat'),1),'g')
        plot(ll_10(SmartVec(1:1:6,10,'flat'),1),'c')
%         pause
    end
    
    % calculate 2-means clustering
    try
        %             if 0
        IDL_1   = kmeans2(ll_1);
        IDL_3   = kmeans2(ll_3);
        IDL_5   = kmeans2(ll_5);
        IDL_10  = kmeans2(ll_10);
        %             else
        %             IDL_0   = kmeans3(ll_0);
        %             IDL_3   = kmeans3(ll_3);
        %             IDL_5   = kmeans3(ll_5);
        %             IDL_10  = kmeans3(ll_10);
        %             end
        
    
    if 0
        % plot original and down-sampled data as well as assignments by
        % kmeans
        ll_1 = bsxfun(@minus,ll_1,mean(ll_1));
        ll_3 = bsxfun(@minus,ll_3,mean(ll_3));
        ll_5 = bsxfun(@minus,ll_5,mean(ll_5));
        ll_10 = bsxfun(@minus,ll_10,mean(ll_10));
        
%         figure(2), clf
        subplot(1,3,2), cla
        hold on
        plot(ll_1(:,1),'b')
        plot(find(IDL_1==IDL_1(1))              ,ll_1(IDL_1==IDL_1(1)),'bo')
        plot(find(IDL_1==3-IDL_1(1))            ,ll_1(IDL_1==3-IDL_1(1)),'ro')
        plot(ll_3(SmartVec(1:1:20,3,'flat'),1),'r')
        plot((find(IDL_3==IDL_3(1))-1)*3+2      ,ll_3(IDL_3==IDL_3(1)),'bo')
        plot((find(IDL_3==3-IDL_3(1))-1)*3+2    ,ll_3(IDL_3==3-IDL_3(1)),'ro')
        plot(ll_5(SmartVec(1:1:12,5,'flat'),1),'g')
        plot((find(IDL_5==IDL_5(1))-1)*5+3      ,ll_5(IDL_5==IDL_5(1)),'bo')
        plot((find(IDL_5==3-IDL_5(1))-1)*5+3    ,ll_5(IDL_5==3-IDL_5(1)),'ro')
        plot(ll_10(SmartVec(1:1:6,10,'flat'),1),'c')
        plot((find(IDL_10==IDL_10(1))-1)*10+5   ,ll_10(IDL_10==IDL_10(1)),'bo')
        plot((find(IDL_10==3-IDL_10(1))-1)*10+5 ,ll_10(IDL_10==3-IDL_10(1)),'ro')
        
%         figure(3), clf
        subplot(1,3,3), cla
        hold on
        plot(ll_1(:,1),ll_1(:,2),'b')
        plot(ll_1(IDL_1==IDL_1(1),1)  ,ll_1(IDL_1==IDL_1(1)  ,2),'bo')
        plot(ll_1(IDL_1==3-IDL_1(1),1),ll_1(IDL_1==3-IDL_1(1),2),'ro')
        plot(ll_3(:,1),ll_3(:,2),'r')
        plot(ll_3(IDL_3==IDL_3(1),1)  ,ll_3(IDL_3==IDL_3(1)  ,2),'bo')
        plot(ll_3(IDL_3==3-IDL_3(1),1),ll_3(IDL_3==3-IDL_3(1),2),'ro')
        plot(ll_5(:,1),ll_5(:,2),'g')
        plot(ll_5(IDL_5==IDL_5(1),1)  ,ll_5(IDL_5==IDL_5(1)  ,2),'bo')
        plot(ll_5(IDL_5==3-IDL_5(1),1),ll_5(IDL_5==3-IDL_5(1),2),'ro')
        plot(ll_10(:,1),ll_10(:,2),'c')
        plot(ll_10(IDL_10==IDL_10(1),1)  ,ll_10(IDL_10==IDL_10(1)  ,2),'bo')
        plot(ll_10(IDL_10==3-IDL_10(1),1),ll_10(IDL_10==3-IDL_10(1),2),'ro')
        
        pause
    end
    catch ER
        if strcmp(ER.identifier,'stats:kmeans:EmptyCluster')
            
            % If an empty cluster error is encountered, try again next
            % iteration. This can occur particularly in long
            % fixations, as the number of clusters there should be 1,
            % but we try to fit 2 to detect a saccade (i.e. 2 fixations)
            
            % visual explanation of empty cluster errors:
            % http://www.ceng.metu.edu.tr/~tcan/ceng465_s1011/Schedule/KMeansEmpty.html
            
            fprintf('Empty cluster error encountered at sample %i. Trying again on next iteration. \n',i);
            counterrors = counterrors + 1;
            continue
        else
            fprintf('Unknown error encountered at sample %i. \n',i);
        end
    end
    
    % detect switches and weight of switch (= 1/number of switches in
    % portion)
    switches1 = abs(diff(IDL_1));
    switches1w = 1/sum(switches1);
    switches3 = abs(diff(IDL_3));
    switches3w = 1/sum(switches3);
    switches5 = abs(diff(IDL_5));
    switches5w = 1/sum(switches5);
    switches10 = abs(diff(IDL_10));
    switches10w = 1/sum(switches10);
    
    % get nearest samples of switch and add weight
    weighted = [switches1*switches1w; 0];
    j = find(switches3)*3;
    for o=-1:1
        weighted(j+o) = weighted(j+o) + switches3w;
    end
    % loop over nswitches instead of window to add to as usually less
    % than 5 (and for lower one definately less than 10) switches in
    % this downsampled data
    for j=find(switches5).'*5
        weighted(j-2:j+2) = weighted(j-2:j+2) + switches5w;
    end
    for j=find(switches10).'*10
        weighted(j-5:j+4) = weighted(j-4:j+5) + switches10w;
    end
    
    % add to totalweights
    totalweights(idx) = totalweights(idx) + weighted;
    % record how many times each sample was tested
    nrtests(idx) = nrtests(idx) + 1;
    
    
    % update i
    i = i + stepsize;
    eind = eind + stepsize;
    qWhichMiss = (on>=i & on<=eind) | (off>=i & off<=eind);
    if any(qWhichMiss)
        % we have some missing in this window. we don't process windows
        % with missing. Move back if we just skipped some samples, or else
        % skip whole missing and place start of window and first next
        % non-missing.
        if on(qWhichMiss)==eind-stepsize+1
            % continue at first non-missing
            i = off(qWhichMiss)+1;
        else
            % we skipped some points, move window back so that we analyze
            % up to first next missing point
            i = min(on(qWhichMiss))-nrsamples;
        end
        eind = i+nrsamples-1;
    end
    if eind>length(xpos) && eind-stepsize<length(xpos)
        % we just exceeded data bound, but previous eind was before end of
        % data: we have some unprocessed samples. retreat just enough so we
        % process those end samples once
        d       = eind-length(xpos);
        eind    = eind-d;
        i       = i-d;
    end
end

% create final weights
finalweights = totalweights./nrtests;

return