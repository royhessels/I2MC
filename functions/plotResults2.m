fprintf('Plotting results \n');

hf = figure(1);
set(hf,'Position',[0 0 1920/2 1080/2])

%% plot layout

myfontsize = 12;
fixlinewidth = 2;

%% plot raw signal

h1 = subplot(3,2,1);
plot(timestamp,llx,'b-');
title('Raw data','FontSize',myfontsize);
ylabel('Horizontal position (pixels)','FontSize',myfontsize);
hold on
plot(timestamp,rrx,'r-');
hold off
axis([0 timestamp(end) 0 xres]);

h3 = subplot(3,2,3);
plot(timestamp,lly,'b-');
ylabel('Vertical position (pixels)','FontSize',myfontsize);
xlabel('Time (ms)','FontSize',myfontsize);
hold on
plot(timestamp,rry,'r-');
hold off
axis([0 timestamp(end) 0 yres]);

%% plot clustering results

h2 = subplot(3,2,2);
plot(timestamp,xpos)
title('Processed data (averaged, interpolated, clustering on average + separate)','FontSize',myfontsize);
axis([0 timestamp(end) 0 xres]);

% add fixations
hold on
for b = 1:length(fixstart)
    plot([fixstart(b) fixend(b)],[xmedian(b) xmedian(b)],'r-','LineWidth',fixlinewidth);
end
hold off

h4 = subplot(3,2,4);
plot(timestamp,ypos)

axis([0 timestamp(end) 0 yres]);

% add fixations
hold on
for b = 1:length(fixstart)
    plot([fixstart(b) fixend(b)],[ymedian(b) ymedian(b)],'r-','LineWidth',fixlinewidth);
end
hold off

h6 = subplot(3,2,6);
plot(timestamp,finalweights)
ylabel('Clustering weight','FontSize',myfontsize);
xlabel('Time (ms)','FontSize',myfontsize);
axis([0 timestamp(end) 0 4]);

% add fixations
hold on
for b = 1:length(fixstart)
    plot([fixstart(b) fixend(b)],[cutoff cutoff],'r-','LineWidth',fixlinewidth);
end
hold off

%% axes stuff Diederick

% link all x axes so zoom (change x limits) on one leads to same new x
% limits for all other ones
linkaxes([h1 h2 h3 h4 h6],'x');
% en nu de truuc, we willen ook de ylims aan elkaar hangen, maar per kolom
% los. Twee keer linkaxes aanroepen werkt niet omdat dan het effect van
% eerdere aanroepen teniet wordt gedaan. Met linkprop kunnen we veel meer..
lnkObj1 = linkprop([h1 h2], 'YLim'); % the link object (return value) 'must exist within the context where you want property linking to occur"
lnkObj2 = linkprop([h3 h4], 'YLim');
lnkObj3 = linkprop([h6], 'YLim');

%% save and close
set(hf,'PaperPositionMode','auto');
drawnow;
% pause;
% print('-depsc2',savefile);
% close(hf);