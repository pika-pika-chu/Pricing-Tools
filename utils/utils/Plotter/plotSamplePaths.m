function plotSamplePaths(idStart,idEnd,paths,baseCcy)


numOfPeriod = width(paths); %#ok<NODEF> 
l = linspace(1,numOfPeriod,numOfPeriod);

figure()
hold on 
plotTitle = ['Sample simulated paths from id #',num2str(idStart),' to id#',num2str(idEnd)];
get(gca,'fontname') ; % shows you what you are using.
set(gca,'fontname','Calibri')  % Set it to times
set(gca,'fontsize',8);
ytickformat(baseCcy);
ytickformat('$%,.0f');
plot(l,paths(idStart:idEnd,:),'LineWidth',1.5)
title(plotTitle)
grid on
hold off



end

