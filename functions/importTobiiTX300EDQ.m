function [timestamp,lx,ly,rx,ry] = importTobiiTX300EDQ(file,res,missingx,missingy)

dat = h5read(file,'/data_collection/events/eyetracker/BinocularEyeSampleEvent');
% "binocular" data consists of same gaze coordinates for both eyes (so just
% take left). Pupil size does seem to be separate for both eyes, but just
% take one anyway
dat = rmfield(dat,setxor(fieldnames(dat),{'time','left_gaze_x','left_gaze_y','right_gaze_x','right_gaze_y','status'}));

msg = h5read(file,'/data_collection/condition_variables/EXP_CV_1');
msg = rmfield(msg,setxor(fieldnames(msg),{'TRIAL_START','TRIAL_END','posx','posy','BLOCK'}));
msg.BLOCK       = num2cell(msg.BLOCK.',2);
msg.BLOCK       = cellfun(@deblank,msg.BLOCK,'uni',false);
qFS = strcmp(msg.BLOCK,'FS');
msg = structfun(@(x) x(qFS),msg,'uni',false);

t1 = msg.TRIAL_START(1);
t2 = msg.TRIAL_END(end);


qTime = dat.time>=t1 & dat.time<=t2;
timestamp   = double(dat.time(qTime)); timestamp = (timestamp-timestamp(1))*1000;
lx          = double(dat.left_gaze_x(qTime))+res(1)/2;
ly          = double(dat.left_gaze_y(qTime))+res(2)/2;
rx          = double(dat.right_gaze_x(qTime))+res(1)/2;
ry          = double(dat.right_gaze_y(qTime))+res(2)/2;
status      = double(dat.status(qTime));

% sometimes we have weird peaks where one sample is (very) far outside the
% monitor. Here, count as missing any data that is more than one monitor
% distance outside the monitor.
qMiss = status==40 | status==44;
lx(qMiss) = missingx;
ly(qMiss) = missingy;
qMiss = status==4 | status==44;
rx(qMiss) = missingx;
ry(qMiss) = missingy;

return