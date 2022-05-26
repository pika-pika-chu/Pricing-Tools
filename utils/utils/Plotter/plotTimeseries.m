function plotTimeseries(ts)
name = ts.getName;
type = char(ts.getType);
type = [upper(type(1)),lower(type(2:end))];
values = ts.getValues;
dates = ts.getDates;
dates = datestr(dates,'mm/dd/yyyy');
dates = datetime(dates);

clf()
leg = [name,' ',type];
get(gca,'fontname') ; % shows you what you are using.
set(gca,'fontname','Calibri')  % Set it to times
set(gca,'fontsize',10);
baseCurerency = lower(name(end-2:end));
ytickformat(baseCurerency);
ytickformat('$%,.0f');




hold on 

plotTitle = [name,' ',type];
plot(dates,values,'LineWidth',1.5)
legend(leg,'Location','northwest')
legend('boxoff')
grid on
title(plotTitle)

x0=10;
y0=10;
width=500;
height=250;
set(gcf,'position',[x0,y0,width,height])
hold off

end

