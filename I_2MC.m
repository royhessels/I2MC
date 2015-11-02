%% IDENTIFICATION BY 2-MEANS CLUSTERING ALGORITHM
% ROY HESSELS - 2014

%% INITIALIZE
clear variables; clear mex; close all; fclose('all'); clc;
dbstop if error;
commandwindow;

%% VARIABLES

% GENERAL
opt.xres                        = 1920; % maximum value of horizontal resolution. used for determing max dispersion and setting plot axes
opt.yres                        = 1080; % maximum value of vertical resolution. used for determing max dispersion and setting plot axes
opt.missingx                    = -opt.xres; % missing value for horizontal position. used throughout functions as signal for data loss
opt.missingy                    = -opt.yres; % missing value for vertical position. used throughout functions as signal for data loss
opt.freq                        = 300; % sampling frequency of data (check that this value matches with values actually obtained from measurement!)

% VISUAL ANGLE STUFF
screendiag                  = 23; % screen diagonal in inches
opt.scrSz                   = [opt.xres opt.yres].*screendiag*2.54./hypot(opt.xres,opt.yres);
% diagpix                     = sqrt(xres^2+yres^2);
% pixpercm                    = diagpix/(screendiag*2.54);
opt.disttoscreen                = 65; % distance to screen in cm. We assume distance to screen is the same; might be slight overestimation though
% rad2deg                     = @(x) x/pi*180;
% degpercm                    = 2*rad2deg(atan(1/(2*disttoscreen)));
% pixperdeg                   = pixpercm/degpercm;

% CUBIC SPLINE INTERPOLATION
opt.windowtimeInterp            = 0.1;  % max duration (s) of missing values for interpolation to occur
opt.edgeSampInterp              = 2;    % amount of data (number of samples) at edges needed for interpolation
opt.maxdisp                     = opt.xres*0.2*sqrt(2); % maximum displacement during missing for interpolation to be possible

% K-MEANS CLUSTERING
opt.windowtime                  = 0.2; % time window (s) over which to calculate 2-means clustering (choose value so that max. 1 saccade can occur)
opt.steptime                    = 0.02;% time window shift (s) for each iteration. Use zero for sample by sample processing
opt.maxerrors                   = 100; % maximum number of errors allowed in k-means clustering procedure before proceeding to next file
opt.downsamples                 = [2 5 10];

% FIXATION DETERMINATION
opt.cutoffstd                   = 2; % number of standard deviations above mean k-means weights will be used as fixation cutoff
opt.maxMergeDist                = 30; % maximum Euclidean distance in pixels between fixations for merging
opt.maxMergeTime                = 30; % maximum time in ms between fixations for merging
opt.minFixDur                   = 40; % minimum fixation duration after merging, fixations with shorter duration are removed from output

%% SET-UP FOLDERS

folders.data                = '../data/TX300'; % folder in which data is stored (each folder in folders.data is considered 1 subject)
folders.output              = 'output'; % folder for output (will use structure in folders.data for saving output)
folders.func                = 'functions'; % folder for functions, will add to matlab path

fs                          = filesep; % platform-independent file seperator

addpath(genpath(folders.func));

%% START ALGORITHM

% create textfile and open for writing fixations
fid = fopen([folders.output fs 'allfixations.txt'],'w');
fprintf(fid,'FixStart\tFixEnd\tFixDur\tXPos\tYPos\tFlankedByDataLoss\tFraction Interpolated\tWeightCutoff\tRMSxy\tBCEA\tFixRangeX\tFixRangeY\tParticipant\tTrial\n');

% go through all folders in 'data' (one folder is assumed to be one
% subject)
[fold,nfold] = FolderFromFolder(folders.data);

for e = 1:nfold     % 3, 16 for one sample peak test
    
    % make output folder
    mkdir(folders.output,fold(e).name);
    
    % go through all files in fold, if fold is empty continue to next
    % folder
    try 
        [file,nfile] = FileFromFolder([folders.data fs fold(e).name],[],'txt');
    catch err  %moet beter: catch specifieke error bij lege folder
        disp('folder is empty, continuing to next folder');
        continue
    end
    
    for f = 1:nfile
        
        %pre-allocate name for saving results to
        savefile = [folders.output fs fold(e).name fs file(f).name(1:end-4)];
        
        %% IMPORT DATA
        
        fprintf('Importing and processing %s/%s \n',fold(e).name,file(f).name)
        
        [timestamp,llx,lly,rrx,rry] = importTobiiTX300([folders.data fs fold(e).name fs file(f).name],1,[opt.xres opt.yres],opt.missingx,opt.missingy);
        
        % check whether we have data, if not, continue to next file
        if isempty(timestamp)
            fprintf('Empty file encountered, continuing to next file \n');
            continue
        end
        
        %% EvENT DETECTION
        data.time = timestamp;
        data.left.X = llx;
        data.left.Y = lly;
        data.right.X = rrx;
        data.right.Y = rry;
        fix = I_2MCfunc(data,opt);
        
        %% PLOT RESULTS
        
        % no function for plotting yet, just a test-specific script
        plotResults2
        
    end
end

%% CLEAN UP

fclose(fid);
