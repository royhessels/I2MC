fprintf('Plotting results \n');

hf = figure(1);
set(hf,'Position',[0 0 1920 1080])

h1 = subplot(3,3,1);
plot(llx,'b-');
title('Horizontal coordinate -  no processing');
hold on
plot(rrx,'r-');
hold off
axis([0 length(llx) 0 xres]);

h4 = subplot(3,3,4);
plot(lly,'b-');
title('Vertical coordinate -  no processing');
hold on
plot(rry,'r-');
hold off
axis([0 length(lly) 0 yres]);

% plot xpos, ypos, finalweights
h2 = subplot(3,3,2);
plot(xpos)
title('Horizontal coordinate -  averaged, interpolated, fixations marked');
axis([0 length(xpos) 0 xres]);

% add fixations
hold on
for b = 1:length(fixstart)
    plot([fixstart(b):1:fixend(b)],repmat(xmedian(b),1,length(fixstart(b):1:fixend(b))),'r-','LineWidth',2);
end
hold off

h5 = subplot(3,3,5);
plot(ypos)
title('Vertical coordinate -  averaged, interpolated, fixations marked');
axis([0 length(ypos) 0 yres]);

% add fixations
hold on
for b = 1:length(fixstart)
    plot([fixstart(b):1:fixend(b)],repmat(ymedian(b),1,length(fixstart(b):1:fixend(b))),'r-','LineWidth',2);
end
hold off

h8 = subplot(3,3,8);
plot(finalweights)
title('2-means clustering weights based on averaged and separate eyes - cutoff for fixations in red');
axis([0 length(finalweights) 0 4]);

% add fixations
hold on
for b = 1:length(fixstart)
    plot([fixstart(b):1:fixend(b)],repmat(cutoff,1,length(fixstart(b):1:fixend(b))),'r-','LineWidth',2);
end
hold off

% ADD SUMMED FINALWEIGHTS TO CHECK WHETHER IT DIFFERS
% plot xpos, ypos, finalweights
h3 = subplot(3,3,3);
plot(xpos)
title('Horizontal coordinate -  averaged, interpolated, fixations marked');
axis([0 length(xpos) 0 xres]);

% add fixations
hold on
for b = 1:length(fixstart2)
    plot([fixstart2(b):1:fixend2(b)],repmat(xmedian2(b),1,length(fixstart2(b):1:fixend2(b))),'r-','LineWidth',2);
end
hold off

h6 = subplot(3,3,6);
plot(ypos)
title('Vertical coordinate -  averaged, interpolated, fixations marked');
axis([0 length(ypos) 0 yres]);

% add fixations
hold on
for b = 1:length(fixstart2)
    plot([fixstart2(b):1:fixend2(b)],repmat(ymedian2(b),1,length(fixstart2(b):1:fixend2(b))),'r-','LineWidth',2);
end
hold off

h9 = subplot(3,3,9);
plot(finalweights_avg)
title('2-means clustering weights based on averaged eyes only - cutoff for fixations in red');
axis([0 length(finalweights_avg) 0 4]);

% add fixations
hold on
for b = 1:length(fixstart2)
    plot([fixstart2(b):1:fixend2(b)],repmat(cutoff2,1,length(fixstart2(b):1:fixend2(b))),'r-','LineWidth',2);
end
hold off

% link all x axes so zoom (change x limits) on one leads to same new x
% limits for all other ones
linkaxes([h1 h2 h3 h4 h5 h6 h8 h9],'x');
% en nu de truuc, we willen ook de ylims aan elkaar hangen, maar per kolom
% los. Twee keer linkaxes aanroepen werkt niet omdat dan het effect van
% eerdere aanroepen teniet wordt gedaan. Met linkprop kunnen we veel meer..
lnkObj1 = linkprop([h1 h2 h3], 'YLim'); % the link object (return value) 'must exist within the context where you want property linking to occur"
lnkObj2 = linkprop([h4 h5 h6], 'YLim');
lnkObj3 = linkprop([h8 h9], 'YLim');

% save and close
set(hf,'PaperPositionMode','auto');
drawnow;
% pause
print('-depsc2',savefile);
close(hf);