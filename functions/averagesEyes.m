function [xpos, ypos, qBMiss, qLMiss, qRMiss, varargout] = averagesEyes(lx,rx,missingx,ly,ry,missingy,varargin)
% Averages data from two eyes. Take one eye if only one was found.

assert(mod(length(varargin),2)==0,'if extra variables given, must be a pair for left and right eye')
varargout = cell(1,length(varargin)/2);

xpos = zeros([length(lx) 1]);
ypos = zeros([length(ly) 1]);
for p=1:2:length(varargin)
    varargout{(p+1)/2} = zeros([length(lx) 1]);
end

% get missing
qLMiss = lx == missingx & ly == missingy;
qRMiss = rx == missingx & ry == missingy;
qBMiss = qLMiss & qRMiss;

q = ~qLMiss & ~qRMiss;
xpos(q) = (lx(q) + rx(q)) ./ 2;
ypos(q) = (ly(q) + ry(q)) ./ 2;
for p=1:2:length(varargin)
    varargout{(p+1)/2} = (varargin{p}(q) + varargin{p+1}(q)) ./ 2;
end

q =  qLMiss & ~qRMiss;
xpos(q) = rx(q);
ypos(q) = ry(q);
for p=1:2:length(varargin)
    varargout{(p+1)/2} = varargin{p+1}(q);
end

q = ~qLMiss & qRMiss;
xpos(q) = lx(q);
ypos(q) = ly(q);
for p=1:2:length(varargin)
    varargout{(p+1)/2} = varargin{p}(q);
end

xpos(qBMiss) = missingx;
ypos(qBMiss) = missingy;
for p=1:2:length(varargin)
    varargout{(p+1)/2} = nan;
end
