clc
clear all
close all

AoM = import_data('./6MWT_AoM_only.csv')
Gyr = import_data('./6MWT_Gyro_only.csv')
Acc = import_data('./6MWT_Acc_only.csv')
AoM_Gyr = import_data('./6MWT_AoM_Gyro.csv')
AoM_Acc = import_data('./6MWT_AoM_Acc.csv')
Gyr_Acc = import_data('./6MWT_Gyro_Acc.csv')
AoM_Gyr_Acc = import_data('./6MWT_AoM_Gyro_Acc.csv')


lw = 1.2
mk = 8
x = [10 20 30 60 90 120 180 240 300 360]
figure(1)
plot(x,AoM.MAE_COT_mean,'r','linewidth',lw)
hold on
plot(x,Gyr.MAE_COT_mean,'m','linewidth',lw)
plot(x,Acc.MAE_COT_mean,'b','linewidth',lw)
% plot(x,AoM_Gyr.MAE_COT_mean,'k','linewidth',lw)
% plot(x,AoM_Acc.MAE_COT_mean,'k','linewidth',lw)
% plot(x,Gyr_Acc.MAE_COT_mean,'k','linewidth',lw)
plot(x,AoM_Gyr_Acc.MAE_COT_mean,'k','linewidth',lw)
h1 = plot(x,AoM.MAE_COT_mean,'rx','linewidth',lw,'markersize',mk)
h2 = plot(x,Gyr.MAE_COT_mean,'mo','linewidth',lw,'markersize',mk)
h3 = plot(x,Acc.MAE_COT_mean,'bd','linewidth',lw,'markersize',mk)
% plot(x,AoM_Gyr.MAE_COT_mean,'k','linewidth',lw,'markersize',mk)
% plot(x,AoM_Acc.MAE_COT_mean,'k','linewidth',lw,'markersize',mk)
% plot(x,Gyr_Acc.MAE_COT_mean,'k','linewidth',lw,'markersize',mk)
h4 = plot(x,AoM_Gyr_Acc.MAE_COT_mean,'k^','linewidth',lw,'markersize',mk)
legend([h1, h2, h3, h4],'AoM','Gyr','Acc','All')
xlabel('Cut off time [sec]')
ylabel('6MWT MAE [m]')
set(gca,'fontsize',15)


figure(2)
plot(x,AoM.MAEP_COT_mean,'r','linewidth',lw)
hold on
plot(x,Gyr.MAEP_COT_mean,'m','linewidth',lw)
plot(x,Acc.MAEP_COT_mean,'b','linewidth',lw)
% plot(x,AoM_Gyr.MAEP_COT_mean,'k','linewidth',lw)
% plot(x,AoM_Acc.MAEP_COT_mean,'k','linewidth',lw)
% plot(x,Gyr_Acc.MAEP_COT_mean,'k','linewidth',lw)
plot(x,AoM_Gyr_Acc.MAEP_COT_mean,'k','linewidth',lw)
h1 = plot(x,AoM.MAEP_COT_mean,'rx','linewidth',lw,'markersize',mk)
h2 = plot(x,Gyr.MAEP_COT_mean,'mo','linewidth',lw,'markersize',mk)
h3 = plot(x,Acc.MAEP_COT_mean,'bd','linewidth',lw,'markersize',mk)
% plot(x,AoM_Gyr.MAEP_COT_mean,'k','linewidth',lw,'markersize',mk)
% plot(x,AoM_Acc.MAEP_COT_mean,'k','linewidth',lw,'markersize',mk)
% plot(x,Gyr_Acc.MAEP_COT_mean,'k','linewidth',lw,'markersize',mk)
h4 = plot(x,AoM_Gyr_Acc.MAEP_COT_mean,'k^','linewidth',lw,'markersize',mk)
legend([h1, h2, h3, h4],'AoM','Gyr','Acc','All')
xlabel('Cut off time [sec]')
ylabel('6MWT MAE% [%]')
set(gca,'fontsize',15)



