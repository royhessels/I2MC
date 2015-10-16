function [finalweights,stopped] = twoClusterWeighting(xpos,ypos,windowtime,freq,missingxvalue,missingyvalue,feedbackiterations,maxerrors)
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
% missingyvalue                 = missing value for vertical position
% feedbackiterations            = number of iterations of MATLAB function kmeans.m are run before giving feedback
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

i=1;
while i<=length(xpos)-nrsamples
    % check if max errors is crossed
    if counterrors > maxerrors
        fprintf('Too many empty clusters encountered, aborting file. \n');
        stopped = 1;
        finalweights = NaN;
        return
    end
    
    if mod(i,feedbackiterations) == 0
        fprintf('Calculating 2-means clustering, now at sample %i \n',i);
    end
    
    weighted = zeros([nrsamples 1]);
    
    % select data portion of nrsamples
    lx = xpos(i:i+(nrsamples-1));
    ly = ypos(i:i+(nrsamples-1));
    
    %%%%% WILL NOW EXCLUDE SAMPLE WITH DATA LOSS
    if sum(lx == missingxvalue) == 0 && sum(ly == missingyvalue) == 0
        
        % downsample to 1/3, 1/5 & 1/10 of original sampling frequency
        
        % note: decimate can also deal with unequal division, i.e. 300/7,
        % unlike averaging over samples
        lx_3 = decimate(lx,3);
        ly_3 = decimate(ly,3);
        lx_5 = decimate(lx,5);
        ly_5 = decimate(ly,5);
        lx_10 = decimate(lx,10);
        ly_10 = decimate(ly,10);
        
        % calculate 2-means clustering
        try
            IDL = kmeans([lx ly],2);
            IDL_3 = kmeans([lx_3 ly_3],2);
            IDL_5 = kmeans([lx_5 ly_5],2);
            IDL_10 = kmeans([lx_10 ly_10],2);
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
        switches0 = abs(diff(IDL));
        switches0w = 1/sum(abs(diff(IDL)));
        switches3 = abs(diff(IDL_3));
        switches3w = 1/sum(abs(diff(IDL_3)));
        switches5 = abs(diff(IDL_5));
        switches5w = 1/sum(abs(diff(IDL_5)));
        switches10 = abs(diff(IDL_10));
        switches10w = 1/sum(abs(diff(IDL_10)));
        
        % get nearest samples of switch and add weight
        for j=find(switches0)
            weighted(j) = weighted(j) + switches0w;
        end
        for j=find(switches3).'
            weighted(j*3:j*3+2) = weighted(j*3:j*3+2) + switches3w;
        end
        for j=find(switches5).'
            weighted(j*5:j*5+4) = weighted(j*5:j*5+4) + switches5w;
        end
        for j=find(switches10).'
            weighted(j*10:j*10+9) = weighted(j*10:j*10+9) + switches10w;
        end
        
        % add to totalweights
        totalweights(i:i+(nrsamples-1)) = totalweights(i:i+(nrsamples-1)) + weighted;
        % record how many times each sample was tested
        nrtests(i:i+(nrsamples-1)) = nrtests(i:i+(nrsamples-1)) + 1;
        
    end
    
    % update i
    i = i + 1;
    
end

% create final weights
finalweights = totalweights./nrtests;

return