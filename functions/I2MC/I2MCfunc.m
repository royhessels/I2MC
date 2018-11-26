function [fix,data,par] = I2MCfunc(data,varargin)
% Hessels, R.S., Niehorster, D.C., Kemner, C., & Hooge, I.T.C., (2016).
% Noise-robust fixation detection in eye-movement data - Identification by 
% 2-means clustering (I2MC). Submitted.

%% deal with inputs
% get inputs the user specified and throw them in the parser
if isstruct(varargin{1})
    % convert to key-value pairs
    assert(isscalar(varargin),'only one input for options is expected if options are given as a struct')
    varargin = [reshape([fieldnames(varargin{1}) struct2cell(varargin{1})].',1,[]) varargin(2:end)];
end

% set defaults
% required parameters:
par.xres            = [];
par.yres            = [];
par.freq            = [];
par.missingx        = [];
par.missingy        = [];
par.scrSz           = [];       % screen size (e.g. in cm). Optional, specify if want fixation statistics in deg
par.disttoscreen    = [];       % screen distance (in same unit as size). Optional, specify if want fixation statistics in deg
% parameters with defaults:
% CUBIC SPLINE INTERPOLATION
par.windowtimeInterp = .1;      % max duration (s) of missing values for interpolation to occur
par.edgeSampInterp  = 2;        % amount of data (number of samples) at edges needed for interpolation
par.maxdisp         = [];       % (default value set below if needed) maximum displacement during missing for interpolation to be possible
% K-MEANS CLUSTERING
par.windowtime      = .2;       % time window (s) over which to calculate 2-means clustering (choose value so that max. 1 saccade can occur)
par.steptime        = .02;      % time window shift (s) for each iteration. Use zero for sample by sample processing
par.downsamples     = [2 5 10]; % downsample levels (can be empty)
par.downsampFilter  = 1;        % use chebychev filter when downsampling? 1: yes, 0: no. requires signal processing toolbox. is what matlab's downsampling functions do, but could cause trouble (ringing) with the hard edges in eye-movement data
par.chebyOrder      = 8;        % order of cheby1 Chebyshev downsampling filter, default is normally ok, as long as there are 25 or more samples in the window (you may have less if your data is of low sampling rate or your window is small
par.maxerrors       = 100;      % maximum number of errors allowed in k-means clustering procedure before proceeding to next file
% FIXATION DETERMINATION
par.cutoffstd       = 2;        % number of standard deviations above mean k-means weights will be used as fixation cutoff
par.onoffsetThresh  = 3;        % number of MAD away from median fixation duration. Will be used to walk forward at fixation starts and backward at fixation ends to refine their placement and stop algorithm from eating into saccades
par.maxMergeDist    = 30;       % maximum Euclidean distance in pixels between fixations for merging
par.maxMergeTime    = 30;       % maximum time in ms between fixations for merging
par.minFixDur       = 40;       % minimum fixation duration (ms) after merging, fixations with shorter duration are removed from output


% loop over input
checkNumeric = @(x,k) assert(isnumeric(x),'The value of ''%s'' is invalid. Expected input to be one of these types:\n\ndouble, single, uint8, uint16, uint32, uint64, int8, int16, int32, int64\n\nInstead its type was %s.',k,class(x));
checkScalar  = @(x,k) assert(isscalar(x),'The value of ''%s'' is invalid. Expected input to be a scalar.',k);
checkNumel2  = @(x,k) assert(numel(x)==2,'The value of ''%s'' is invalid. Expected input to be an array with number of elements equal to 2.',k);
checkInt     = @(x,k) assert(~any(mod(x,1)),'The value of ''%s'' is invalid. Expected input to be integer-valued.',k);
for p=1:2:length(varargin)
    key = varargin{p};
    if p+1>length(varargin)
        error('No value was given for ''%s''. Name-value pair arguments require a name followed by a value.',key);
    end
    value = varargin{p+1};
    switch key
        case {'xres','yres','freq','missingx','missingy','disttoscreen','windowtimeInterp','maxdisp','windowtime','steptime','cutoffstd','onoffsetThresh','maxMergeDist','maxMergeTime','minFixDur'}
            checkNumeric(value,key);
            checkScalar(value,key);
            par.(key) = value;
        case {'downsampFilter','chebyOrder','maxerrors','edgeSampInterp'}
            checkInt(value,key);
            checkScalar(value,key);
            par.(key) = value;
        case 'scrSz'
            checkNumeric(value,key);
            checkNumel2(value,key);
            par.(key) = value;
        case 'downsamples'
            checkInt(value,key);
            par.(key) = value;
        otherwise
            if ~ischar(key)
                error('Expected a string for the parameter name at position %d, instead the input type was ''%s''.',class(key));
            else
                error('Key "%s" not recognized',key);
            end
    end
end

% deal with required options
% if empty, user did not specify these
checkFun = @(opt,str) assert(~isempty(par.(opt)),'I2MCfunc: %s must be specified using the ''%s'' option',str,opt);
checkFun('xres', 'horizontal screen resolution')
checkFun('yres',   'vertical screen resolution')
checkFun('freq', 'tracker sampling rate')
checkFun('missingx', 'value indicating data loss for horizontal position')
checkFun('missingy', 'value indicating data loss for vertical position')
% process parameters with defaults
if isempty(par.maxdisp)
    par.maxdisp               = par.xres*0.2*sqrt(2); % maximum displacement during missing for interpolation to be possible
end

% check filter
if par.downsampFilter
    assert(exist('cheby1','file')==2,'I2MCfunc: When setting the ''downsampFilter'' option to true, the function ''cheby1'' from the signal processing toolbox is required. It appears this function is not available in your installation. Set the option to 0.')
    nSampRequired = max(1,3*par.chebyOrder)+1;  % nSampRequired = max(1,3*(nfilt-1))+1, where nfilt = chebyOrder+1
    nSampInWin    = round(par.windowtime/(1/par.freq));
    assert(nSampInWin>=nSampRequired,'I2MCfunc: Filter parameters requested with the setting ''chebyOrder'' will not work for the sampling frequency of your data. Please lower ''chebyOrder'', or set the setting ''downsampFilter'' to 0')
end
assert(~any(mod(par.freq,par.downsamples)),'I2MCfunc: Some of your downsample levels are not divisors of your sampling frequency. Change the option ''downsamples''')

% setup visual angle conversion
pixperdeg = [];
if ~isempty(par.scrSz) && ~isempty(par.disttoscreen)
    pixpercm                    = mean([par.xres par.yres]./par.scrSz(:).');
    rad2deg                     = @(x) x/pi*180;
    degpercm                    = 2*rad2deg(atan(1/(2*par.disttoscreen)));
    pixperdeg                   = pixpercm/degpercm;
end

%% START ALGORITHM

%% PREPARE INPUT DATA
% make sure all fields in data are columns
fs = fieldnames(data);
for f=1:length(fs)
    if isstruct(data.(fs{f}))
        fs2 = fieldnames(data.(fs{f}));
        for f2=1:length(fs2)
            data.(fs{f}).(fs2{f2}) = data.(fs{f}).(fs2{f2})(:);
        end
    else
        data.(fs{f}) = data.(fs{f})(:);
    end
end
% deal with monocular data, or create average over two eyes
if isfield(data,'left') && ~isfield(data,'right')
    xpos = data.left.X;
    ypos = data.left.Y;
    missing = isnan(data.left.X) | data.left.X==par.missingx | isnan(data.left.Y) | data.left.Y==par.missingy;
    data.left.missing = missing;
    q2Eyes = false;
elseif isfield(data,'right') && ~isfield(data,'left')
    xpos = data.right.X;
    ypos = data.right.Y;
    missing = isnan(data.right.X) | data.right.X==par.missingx | isnan(data.right.Y) | data.right.Y==par.missingy;
    data.right.missing = missing;
    q2Eyes = false;
elseif isfield(data,'average')
    xpos = data.average.X;
    ypos = data.average.Y;
    missing = isnan(data.average.X) | data.average.X==par.missingx | isnan(data.average.Y) | data.average.Y==par.missingy;
    data.average.missing = missing;
    q2Eyes = isfield(data,'right') && isfield(data,'left');
    if q2Eyes
        % we have left and right and average already provided, but we need
        % to get missing in the individual eye signals
        [llmiss, rrmiss] = getMissing(data.left.X,data.right.X,par.missingx,data.left.Y,data.right.Y,par.missingy);
        data.left.missing  = llmiss;
        data.right.missing = rrmiss;
    end
else % we have left and right, average them
    [data.average.X, data.average.Y, missing, llmiss, rrmiss] = averageEyes(data.left.X,data.right.X,par.missingx,data.left.Y,data.right.Y,par.missingy);
    xpos = data.average.X;
    ypos = data.average.Y;
    data.average.missing = missing;
    data.left.missing    = llmiss;
    data.right.missing   = rrmiss;
    q2Eyes = true;
end

%% INTERPOLATION

% get interpolation windows for average and individual eye signals
fprintf('Searching for valid interpolation windows\n');
interpwins   = findInterpWins(xpos, ypos, missing,par.windowtimeInterp,par.edgeSampInterp,par.freq,par.maxdisp);
if q2Eyes
    llinterpwins = findInterpWins(data.left.X  ,data.left.Y  ,llmiss ,par.windowtimeInterp,par.edgeSampInterp,par.freq,par.maxdisp);
    rrinterpwins = findInterpWins(data.right.X ,data.right.Y ,rrmiss ,par.windowtimeInterp,par.edgeSampInterp,par.freq,par.maxdisp);
end

% Use Steffen interpolation and replace values
fprintf('Replace interpolation windows with Steffen interpolation\n');
[xpos,ypos,missingn]= windowedInterpolate(xpos, ypos, missing, interpwins,par.edgeSampInterp);
if q2Eyes
    [llx ,lly ,llmiss]  = windowedInterpolate(data.left.X  ,data.left.Y  ,llmiss ,llinterpwins,par.edgeSampInterp);
    [rrx ,rry ,rrmiss]  = windowedInterpolate(data.right.X ,data.right.Y ,rrmiss ,rrinterpwins,par.edgeSampInterp);
end


if ~q2Eyes
    %% CALCULATE 2-MEANS CLUSTERING FOR SINGLE EYE
    
    % get kmeans-clustering for averaged signal
    fprintf('2-Means clustering started for averaged signal \n');
    [data.finalweights,stopped] = twoClusterWeighting(xpos,ypos,missingn,par.downsamples,par.downsampFilter,par.chebyOrder,par.windowtime,par.steptime,par.freq,par.maxerrors);
    
    % check whether clustering succeeded
    if stopped
        fprintf('Clustering stopped after exceeding max errors, continuing to next file \n');
        return
    end
    
    %% CALCULATE 2-MEANS CLUSTERING FOR SEPARATE EYES
elseif q2Eyes
    % get kmeans-clustering for left eye signal
    fprintf('2-Means clustering started for left eye signal \n');
    [finalweights_left,stopped] = twoClusterWeighting(llx,lly,llmiss,par.downsamples,par.downsampFilter,par.chebyOrder,par.windowtime,par.steptime,par.freq,par.maxerrors);
    
    % check whether clustering succeeded
    if stopped
        fprintf('Clustering stopped after exceeding max errors, continuing to next file \n');
        return
    end
    
    % get kmeans-clustering for right eye signal
    fprintf('2-Means clustering started for right eye signal \n');
    [finalweights_right,stopped] = twoClusterWeighting(rrx,rry,rrmiss,par.downsamples,par.downsampFilter,par.chebyOrder,par.windowtime,par.steptime,par.freq,par.maxerrors);
    
    % check whether clustering succeeded
    if stopped
        fprintf('Clustering stopped after exceeding max errors, continuing to next file \n');
        return
    end
    
    %% AVERAGE FINALWEIGHTS OVER COMBINED & SEPARATE EYES
    data.finalweights = mean([finalweights_left finalweights_right],2,'omitnan');
end

%% DETERMINE FIXATIONS BASED ON FINALWEIGHTS_AVG
fprintf('Determining fixations based on clustering weight mean for averaged signal and separate eyes + 2*std \n')
[fix.cutoff,fix.start,fix.end,fix.startT,fix.endT,fix.dur,fix.xpos,fix.ypos,fix.flankdataloss,fix.fracinterped] = getFixations(data.finalweights,data.time,xpos,ypos,missing,par.cutoffstd,par.onoffsetThresh,par.maxMergeDist,par.maxMergeTime,par.minFixDur);
[fix.RMSxy,fix.BCEA,fix.fixRangeX,fix.fixRangeY] = getFixStats(xpos,ypos,missing,fix.start,fix.end,pixperdeg);
