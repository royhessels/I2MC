function [finalweights,stopped] = twoClusterWeightingNoLossOpt(xpos,ypos,windowtime,freq,missingxvalue,missingyvalue,maxerrors)
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
% freq                          = sampling frequency of data
% missingxvalue                 = missing value for horizontal position
% missingxvalue                 = missing value for horizontal position
% maxerrors                     = maximum number of errors allowed in k-means clustering procedure before proceeding to next file

% Output:
% finalweights                  = vector of 2-means clustering weights (one weight for each sample), the higher, the more likely a saccade happened
% stopped                       = whether maxerrors was reached or not

% Roy Hessels - 2014

% calculate number of samples of the moving window
nrsamples = round(windowtime/(1/freq));

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
idx_3 = fliplr(nrsamples:-3:1);
idx_5 = fliplr(nrsamples:-5:1);
idx_10 = fliplr(nrsamples:-10:1);

% see where are missing in this data, for better running over the data
% below
[on,off] = bool2bounds(xpos == missingxvalue | ypos == missingyvalue);
if ~isempty(on)&&on(1)==1
    i=off(1)+1;  % start at first non-missing
else
    i=1;
end
eind = i+nrsamples-1;
while i<=length(xpos)-nrsamples
    % check if max errors is crossed
    if counterrors > maxerrors
        fprintf('Too many empty clusters encountered, aborting file. \n');
        stopped = 1;
        finalweights = NaN;
        return
    end
    
    % select data portion of nrsamples
    idx = i:eind;
    % the filtering done in decimate is done once above. now we simply need
    % to select each nth sample where n is the integer factor by which
    % number of samples is reduced. select samples such that they are till
    % end of window
    
    %%%%% WILL NOW EXCLUDE SAMPLE WITH DATA LOSS
    ll_0 = [xpos(idx) ypos(idx)];
    qWhichMiss = (off>i & off<eind) | (on>i & on<eind);
    if any(qWhichMiss)
        % continue at first non-missing
        i = max(off(qWhichMiss))+1;
        eind = i+nrsamples-1;
        continue;
    end
    
    ll_3 = filtfilt(b3,a3,ll_0);
    ll_3 = ll_3(idx_3,:);
    ll_5 = filtfilt(b5,a5,ll_0);
    ll_5 = ll_5(idx_5,:);
    ll_10 = filtfilt(b10,a10,ll_0);
    ll_10 = ll_10(idx_10,:);
    
    % calculate 2-means clustering
    try
        %             if 0
        IDL_0   = kmeans2(ll_0);
        IDL_3   = kmeans2(ll_3);
        IDL_5   = kmeans2(ll_5);
        IDL_10  = kmeans2(ll_10);
        %             else
        %             IDL_0   = kmeans3(ll_0);
        %             IDL_3   = kmeans3(ll_3);
        %             IDL_5   = kmeans3(ll_5);
        %             IDL_10  = kmeans3(ll_10);
        %             end
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
    switches0 = abs(diff(IDL_0));
    switches0w = 1/sum(switches0);
    switches3 = abs(diff(IDL_3));
    switches3w = 1/sum(switches3);
    switches5 = abs(diff(IDL_5));
    switches5w = 1/sum(switches5);
    switches10 = abs(diff(IDL_10));
    switches10w = 1/sum(switches10);
    
    % get nearest samples of switch and add weight
    weighted = [switches0*switches0w; 0];
    j = find(switches3)*3;
    for o=0:2
        weighted(j+o) = weighted(j+o) + switches3w;
    end
    % loop over nswitches instead of window to add to as usually less
    % than 5 (and for lower one definately less than 10) switches in
    % this downsampled data
    for j=find(switches5).'*5
        weighted(j:j+4) = weighted(j:j+4) + switches5w;
    end
    for j=find(switches10).'*10
        weighted(j:j+9) = weighted(j:j+9) + switches10w;
    end
    
    % add to totalweights
    totalweights(idx) = totalweights(idx) + weighted;
    % record how many times each sample was tested
    nrtests(idx) = nrtests(idx) + 1;
    
    
    % update i
    i = i + 1;
    eind = eind + 1;
end

% create final weights
finalweights = totalweights./nrtests;

return