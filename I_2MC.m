%% IDENTIFICATION BY 2-MEANS CLUSTERING ALGORITHM
% ROY HESSELS - 2014

%% INITIALIZE
clear variables; clear mex; close all; fclose('all'); clc;
dbstop if error;
commandwindow;

%% VARIABLES

% GENERAL
xres                        = 1920; % maximum value of horizontal resolution. used for determing max dispersion and setting plot axes
yres                        = 1080; % maximum value of vertical resolution. used for determing max dispersion and setting plot axes
missingx                    = -xres; % missing value for horizontal position. used throughout functions as signal for data loss
missingy                    = -yres; % missing value for vertical position. used throughout functions as signal for data loss
freq                        = 300; % sampling frequency of data (check that this value matches with values actually obtained from measurement!)

% VISUAL ANGLE STUFF
screendiag                  = 23; % screen diagonal in inches
diagpix                     = sqrt(xres^2+yres^2);
pixpercm                    = diagpix/(screendiag*2.54);
disttoscreen                = 65; % distance to screen in cm. We assume distance to screen is the same; might be slight overestimation though
rad2deg                     = @(x) x/pi*180;
degpercm                    = 2*rad2deg(atan(1/(2*disttoscreen)));
pixperdeg                   = pixpercm/degpercm;

% CUBIC SPLINE INTERPOLATION
windowtimeInterp            = 0.1;  % max duration (s) of missing values for interpolation to occur
edgeSampInterp              = 2;    % amount of data (number of samples) at edges needed for interpolation
maxdisp                     = xres*0.2*sqrt(2); % maximum displacement during missing for interpolation to be possible

% K-MEANS CLUSTERING
windowtime                  = 0.2; % time window (s) over which to calculate 2-means clustering (choose value so that max. 1 saccade can occur)
steptime                    = 0.02;% time window shift (s) for each iteration. Use zero for sample by sample processing
maxerrors                   = 100; % maximum number of errors allowed in k-means clustering procedure before proceeding to next file

% FIXATION DETERMINATION
cutoffstd                   = 2; % number of standard deviations above mean k-means weights will be used as fixation cutoff
maxMergeDist                = 30; % maximum Euclidean distance in pixels between fixations for merging
maxMergeTime                = 30; % maximum time in ms between fixations for merging
minFixDur                   = 40; % minimum fixation duration after merging, fixations with shorter duration are removed from output

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

for e = 1%:nfold     % 3, 16 for one sample peak test
    
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
    
    for f = 8%1:nfile
        
        %pre-allocate name for saving results to
        savefile = [folders.output fs fold(e).name fs file(f).name(1:end-4)];
        
        %% IMPORT DATA
        
        fprintf('Importing and processing %s/%s \n',fold(e).name,file(f).name)
        
        [timestamp,llx,lly,rrx,rry] = importTobiiTX300([folders.data fs fold(e).name fs file(f).name],1,[xres yres],missingx,missingy);
        
        % check whether we have data, if not, continue to next file
        if isempty(timestamp)
            fprintf('Empty file encountered, continuing to next file \n');
            continue
        end
        
        %% CREATE AVERAGE OVER TWO EYES
        [xpos, ypos, missing, llmiss, rrmiss] = averagesEyes(llx,rrx,missingx,lly,rry,missingy);
        
        %% CUBIC SPLINE INTERPOLATION
        % Frank, Vul, & Johnson (2009) interpolate over 100 ms using cubic spline.
        % Saez de Urabain, Johnson, & Smith (2014) interpolate only in fixations
        % based on given time and dispersion between adjacent fixations.
        % We interpolate using Steffen interpolation.
        
        % get interpolation windows for average and individual eye signals
        fprintf('Searching for valid interpolation windows\n');
        interpwins   = findInterpWins(xpos,ypos,missing,windowtimeInterp,edgeSampInterp,freq,maxdisp);
        llinterpwins = findInterpWins(llx ,lly ,llmiss ,windowtimeInterp,edgeSampInterp,freq,maxdisp);
        rrinterpwins = findInterpWins(rrx ,rry ,rrmiss ,windowtimeInterp,edgeSampInterp,freq,maxdisp);
        
        % Use Steffen interpolation and replace values
        fprintf('Replace interpolation windows with Steffen interpolation\n');
        [xpos,ypos,missingn]= windowedInterpolate(xpos,ypos,missing,  interpwins,edgeSampInterp);
        [llxI,llyI]  = windowedInterpolate(llx ,lly ,llmiss ,llinterpwins,edgeSampInterp);
        [rrxI,rryI]  = windowedInterpolate(rrx ,rry ,rrmiss ,rrinterpwins,edgeSampInterp);
        [llxS,llyS,llmiss]  = windowedInterpolate(llx ,lly ,llmiss ,llinterpwins,edgeSampInterp,true);
        [rrxS,rryS,rrmiss]  = windowedInterpolate(rrx ,rry ,rrmiss ,rrinterpwins,edgeSampInterp,true);
        figure(2), hold on;
        plot(timestamp,llx)
        plot(timestamp,llxS,'g')
        plot(timestamp,llxI,'r')
        title('X')
        figure(3), hold on;
        plot(timestamp,lly)
        plot(timestamp,llyS,'g')
        plot(timestamp,llyI,'r')
        title('Y')
        
        %% CALCULATE 2-MEANS CLUSTERING FOR AVERAGED EYES
        % get kmeans-clustering for averaged signal
        fprintf('2-Means clustering started for averaged signal \n');
        % version optimized further, walking over data in bigger steps
        [finalweights,stopped] = twoClusterWeighting2(xpos,ypos,missingn,windowtime,steptime,freq,maxerrors);
        % optimized version of the original, but producing exactly same
        % result (if same random seed, see help rng) as original
%         [finalweights2,stopped] = twoClusterWeightingNoLossOpt(xpos,ypos,windowtime,freq,missingx,missingy,maxerrors);
        % should be pretty much equivalent to:
%         [finalweights2,stopped] = twoClusterWeighting2(xpos,ypos,windowtime,0,freq,missingx,missingy,maxerrors);
        % finally, this version has a different location of the downsampled
        % value (choose from middle of window instead of end), as well as
        % different way of placing switches into weight vector (centered on
        % the downsampled data point). It makes no real difference at all
%         [finalweights,stopped] = twoClusterWeighting3(xpos,ypos,windowtime,steptime,freq,missingx,missingy,maxerrors);
        
%         clf
%         subplot(3,1,1)
%         plot(1:length(xpos),xpos)
%         subplot(3,1,2)
%         plot(1:length(xpos),ypos)
%         subplot(3,1,3)
%         plot(1:length(xpos),finalweights)
%         hold on
%         plot(1:length(xpos),finalweights2,'r')
%         linkaxes([subplot(3,1,1) subplot(3,1,2) subplot(3,1,3)],'x')
%         subplot(3,1,1), axis tight
%         pause
        
        % check whether clustering succeeded
        if stopped
            fprintf('Clustering stopped after exceeding max errors, continuing to next file \n');
            continue
        end
        
        %% CALCULATE 2-MEANS CLUSTERING FOR SEPARATE EYES
        
        % get kmeans-clustering for left eye signal
        fprintf('2-Means clustering started for left eye signal \n');
        [finalweights_left,stopped] = twoClusterWeighting2(llxI,llyI,llmiss,windowtime,steptime,freq,maxerrors);
        
        % check whether clustering succeeded
        if stopped
            fprintf('Clustering stopped after exceeding max errors, continuing to next file \n');
            continue
        end
        
        % get kmeans-clustering for right eye signal
        fprintf('2-Means clustering started for right eye signal \n');
        [finalweights_right,stopped] = twoClusterWeighting2(rrxI,rryI,rrmiss,windowtime,steptime,freq,maxerrors);
        
        % check whether clustering succeeded
        if stopped
            fprintf('Clustering stopped after exceeding max errors, continuing to next file \n');
            continue
        end
        
        %% AVERAGE FINALWEIGHTS OVER COMBINED & SEPARATE EYES
        
        finalweights_avg = nanmean([finalweights finalweights_left finalweights_right],2);
        
        %% DETERMINE FIXATIONS BASED ON FINALWEIGHTS - HAS NOW BECOME OBSOLETE?
        
%         fprintf('Determining fixations based on clustering weight mean for averaged signal + 2*std \n')
%         
%         [cutoff2,fixstart2,fixend2,fixdur2,xmedian2,ymedian2,flankdataloss2,~,~,~,~] = getFixations(finalweights,cutoffstd,xpos,ypos,timestamp,pixperdeg,maxMergeDist,maxMergeTime);
        
        %% DETERMINE FIXATIONS BASED ON FINALWEIGHTS_AVG
        
        fprintf('Determining fixations based on clustering weight mean for averaged signal and separate eyes + 2*std \n')
        
        [cutoff,fixstart,fixend,fixdur,xmedian,ymedian,flankdataloss,fracinterped,RMSxy,BCEA,fixRangeX,fixRangeY] = getFixations(finalweights_avg,missing,xpos,ypos,timestamp,cutoffstd,pixperdeg,maxMergeDist,maxMergeTime,minFixDur);
        
        %% SAVE FIXATIONS (BASED ON COMBINED + SEPARATE EYES) TO FILE
        
        for g=1:numel(fixstart)
            fprintf(fid,'%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\t%s\n',[fixstart(g) fixend(g) fixdur(g) xmedian(g) ymedian(g) flankdataloss(g) fracinterped(g) cutoff, RMSxy(g), BCEA(g), fixRangeX(g), fixRangeY(g)],fold(e).name,file(f).name(1:end-4));
        end
        
        %% PLOT RESULTS
        
        % no function for plotting yet, just a test-specific script
        plotResults2
        
    end
end

%% CLEAN UP

fclose(fid);
