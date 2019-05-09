function [] = nwa_plot_ppi(PHYS1,PHYS2,PSY,varargin)

%%  plot info
nr = 3;
nc = 2;
npplot = nr*nc;
sfac = 6;

%% prep the input

% demean time series
PHYS1 = PHYS1 - mean(PHYS1);
PHYS2 = PHYS2 - mean(PHYS2);

% decenter task regressor
task1 = PSY>=0.1;
task2 = PSY<0.1;
maxPSY      = max(PSY);
minPSY      = min(PSY);
subval      = mean([maxPSY;minPSY]);
PSY         = PSY - subval;


% calc PPI
if length(varargin)<1
    PPI1 = PHYS1.*PSY;
    PPI2 = PHYS2.*PSY;
end

% interaction between times series
PHYSI = time_varying_estimate('tukey',([PHYS1 PHYS2]),sfac);

%% plot
figure('Name','PPI plot')

% 
subplot(nr,nc,1)
plot(PHYS1,'r'); hold on;
plot(PHYS2,'b')
plot(PSY,'k');
plot(PHYSI,'g')
legend('phys1','phys2','PSY','PHYSI')
xlabel('Time')
title('Time-series & Task')

subplot(nr,nc,2)
plot(zscore(PHYS1),'r'); hold on;
plot(zscore(PHYS2),'b')
plot(zscore(PPI1),'k');
plot(zscore(PPI2),'g');
legend('phys1','phys2','PPI1','PPI2')
xlabel('Time')
title('Time-series & PPI terms')

subplot(nr,nc,3)
scatter(PHYS1(task1),PHYS2(task1),'g')
hold on
scatter(PHYS1(task2),PHYS2(task2),'k')

legend('task1','task2')
xlabel('PHYS1')
ylabel('PHYS2')

[coef_fit1,s] = polyfit(PHYS1(task1),PHYS2(task1),1);
[yfit1,dy] = polyconf(coef_fit1,xlim,s,'predopt','curve');
plot(xlim,yfit1,'color','g','LineWidth',2)

[coef_fit2,s] = polyfit(PHYS1(task2),PHYS2(task2),1);
[yfit2,dy] = polyconf(coef_fit2,xlim,s,'predopt','curve');
hold on
plot(xlim,yfit2,'color','k','LineWidth',2)

title('Time-series scatter')

subplot(nr,nc,4)
dat11 = PHYS1(task1);
dat12 = PHYS1(task2);
dat21 = PHYS2(task1);
dat22 = PHYS2(task2);
minphys = min([PHYS1;PHYS2]);
maxphys = max([PHYS1;PHYS2]);
steps = ((maxphys-minphys)/(sum(task2)))*5;
b = minphys:steps:maxphys;
[h11 x11] = histc(dat11,b);
[h12 x12] = histc(dat12,b);
[h21 x21] = histc(dat21,b);
[h22 x22] = histc(dat22,b);
p = plot(b,h11,'r');
hold on
p = plot(b,h12,'--r');
p = plot(b,h21,'b');
p = plot(b,h22,'--b');
title('Time-series distribution')
legend('PHYS1 task1','PHYS1 task2','PHYS2 task1','PHYS2 task2')

subplot(nr,nc,5)
C = corr([PHYS1 PHYS2 PPI1 PPI2 PSY PHYSI]);
nwa_plot_conn(C,'varnames',{'PHYS1' 'PHYS2' 'PPI1' 'PPI2' 'PSY' 'PHYSI'},'title','First-order Correlations')

subplot(nr,nc,6)
Cp = partialcorr([PHYS1 PHYS2 PPI1 PPI2 PSY PHYSI]);
nwa_plot_conn(Cp,'varnames',{'PHYS1' 'PHYS2' 'PPI1' 'PPI2' 'PSY' 'PHYSI'},'title','Partial Correlations')


figure('Name','PPI')
lw = 2;
subplot(6,1,1)
plot(PSY,'k','LineWidth',lw); %ylabel('Task')
subplot(6,1,2)
plot(PHYS1,'r','LineWidth',lw); %ylabel('Time Series R1')
% hYLabel = get(gca,'YLabel');
% %set(hYLabel,'rotation',0,'VerticalAlignment','middle')
% 
% ylp = get(hYLabel,'Position');
% ext=get(hYLabel,'Extent');
% set(hYLabel, 'Rotation',0, 'Position',ylp+[ext(3) 0 0])
subplot(6,1,3)
plot(PPI1,'b','LineWidth',lw); %ylabel('PPI')
% calc the resiudals 
subplot(6,1,4)
stats = regstats(PHYS1,[PSY PPI1]);
plot(stats.r,'g','LineWidth',lw); %ylabel('R1 Residuals')

subplot(6,1,5)
plot(PHYS2,'-.r','LineWidth',lw); %ylabel('Time Series N2')

subplot(6,1,6)
stats = regstats(PHYS2,[PSY PPI2]);
plot(stats.r,'-.g','LineWidth',lw); %ylabel('R2 Residuals')


end

