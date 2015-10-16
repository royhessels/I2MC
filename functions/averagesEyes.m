function [xpos, ypos, qBMiss, qLMiss, qRMiss] = averagesEyes(lx,rx,missingx,ly,ry,missingy)
% Averages data from two eyes. Take one eye if only one was found.

xpos = zeros([length(lx) 1]);
ypos = zeros([length(ly) 1]);

if 0
    xpos(lx ~= missingx & rx ~= missingx) = (lx(lx ~= missingx & rx ~= missingx) + rx(lx ~= missingx & rx ~= missingx)) ./ 2;
    xpos(lx == missingx & rx ~= missingx) = rx(lx == missingx & rx ~= missingx);
    xpos(lx ~= missingx & rx == missingx) = lx(lx ~= missingx & rx == missingx);
    xpos(lx == missingx & rx == missingx) = missingx;
    ypos(ly ~= missingy & ry ~= missingy) = (ly(ly ~= missingy & ry ~= missingy) + ry(ly ~= missingy & ry ~= missingy)) ./ 2;
    ypos(ly == missingy & ry ~= missingy) = ry(ly == missingy & ry ~= missingy);
    ypos(ly ~= missingy & ry == missingy) = ly(ly ~= missingy & ry == missingy);
    ypos(ly == missingy & ry == missingy) = missingy;
else
    % if x is missing, y is also always missing (and corolary also true)
    qLMiss = lx == missingx;
    qRMiss = rx == missingx;
    qBMiss = qLMiss & qRMiss;
    
    q = ~qLMiss & ~qRMiss;
    xpos(q) = (lx(q) + rx(q)) ./ 2;
    ypos(q) = (ly(q) + ry(q)) ./ 2;
    
    q =  qLMiss & ~qRMiss;
    xpos(q) = rx(q);
    ypos(q) = ry(q);
    
    q = ~qLMiss & qRMiss;
    xpos(q) = lx(q);
    ypos(q) = ly(q);
    
    xpos(qBMiss) = missingx;
    ypos(qBMiss) = missingy;
end