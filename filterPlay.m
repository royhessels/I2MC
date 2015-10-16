r = 3;
n = 16:36500;
lx = length(n);
x = sin(2*pi*n/1530) + cos(2*pi*n/1270);

plot(0:lx-1,x)
hold on

x1 = downsample(x,r);
x2 = decimate(x,r);
x3 = decimate(x,r,13);
plot(0:r:lx-1,x1,'ok')
plot(lx-1:-r:0,fliplr(x2),'or')
plot(lx-1:-r:0,fliplr(x3),'og')

% generate pink noise:
pnG = dsp.ColoredNoise(1,length(n),1);
pn = step(pnG).'/5;

figure,plot(pn)

figure
plot(0:lx-1,x,'k')
hold on
xn = x+pn;
plot(0:lx-1,xn,'b')

x2 = decimate(xn,r);
plot(lx-1:-r:0,fliplr(x2),'or-')

% do only the filter, not the decimate
filtmag_db = @(b,a,f) 20*log10(abs((exp(-1i*(0:length(b)-1)*pi*f)*b(:))/(exp(-1i*(0:length(a)-1)*pi*f)*a(:))));
rip = .05;	% passband ripple in dB

[b,a] = cheby1(8, rip, .8/r);
while all(b==0) || (abs(filtmag_db(b,a,.8/r)+rip)>1e-6)
    nfilt = nfilt - 1
    if nfilt == 0
        break
    end
    [b,a] = cheby1(nfilt, rip, .8/r);
end
y1 = filtfilt(b,a,xn);
plot(0:lx-1,y1,'g')

if 0
    % conclusion for below: for the filters we're designing this makes no
    % practical difference, simply use b,a form and filtfilt
    % try fex filtfilthd that takes a dfilt object, as b,a representation can
    % suffer from numerical precision errors. If no real difference in speed,
    % use this better way
    % see first answer here:
    % http://stackoverflow.com/questions/5591278/high-pass-filtering-in-matlab
    [z,p,k]=cheby1(8, rip, .8/r);
    [s,g]=zp2sos(z,p,k);%# create second order sections
    Hd=dfilt.df2sos(s,g);
    y2 = filtfilthd(Hd,xn);
    
    plot(0:lx-1,y2,'c')
end