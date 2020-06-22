clc
clear all
close all


LR = [5	60.6	16.47	48.99	40.73;
10	52.47	13.42	42.09	33.7;
20	48.09	12.21	37.5	30.69;
30	48.73	12.37	37.5	30.37;
60	45.99	10.53	36.24	27.49;
90	44.03	9.87	34.87	25.88;
120	43.05	8.89	34.65	24.69;
180	41.72	8.71	33.36	24.49;
240	40.07	8.39	31.93	23.75;
300	38.61	8.65	30.57	23.6;
360	37.55	8.29	29.43	23.56]

SVR_RBF = [5	60.19	15.04	48.65	37.38;
10	53.09	12.35	43.63	31.52;
20	52.76	11.76	44.04	30.2;
30	52.46	12.47	43.64	29.71;
60	49.82	13.41	40.73	27.93;
90	49.58	11.81	42.09	26.21;
120	49.22	11.3	41.73	25.67;
180	44.51	11.97	35.42	26.56;
240	44.75	12.02	36.51	25.28;
300	48.47	10.9	39.92	25.4;
360	46.97	11.43	37.81	25.96]

SVR_LR = [5	61.09	15	49.47	37.62;
10	54.01	14.24	43.31	33.65;
20	50.04	13.22	40.53	31.34;
30	48.07	13.25	38.25	31.57;
60	46.8	14.2	36.33	30.03;
90	44.38	13.53	34.09	28.74;
120	45.49	12.76	34.5	29.63;
180	45.12	12.81	33.89	29.28;
240	44.75	10.45	33.66	28.75;
300	43.2	11.7	30.8	28.87;
360	42.07	11.89	29.71	28.91]

lw = 2
mk = 10;

x = LR(:,1)
y = LR(:,2)
err = LR(:,3)
figure(1)
errorbar(x,y,err,'color',[0 0 0],'linewidth',2)
hold on

x =SVR_RBF(:,1)
y = SVR_RBF(:,2)
err = SVR_RBF(:,3)
figure(1)
errorbar(x,y,err,'color',[0 0 1],'linewidth',2)

x =SVR_LR(:,1)
y = SVR_LR(:,2)
err = SVR_LR(:,3)
figure(1)
errorbar(x,y,err,'color',[1 0 0],'linewidth',2)
xlabel('Time [Sec]')
ylabel('RMSE [m]')
set(gca,'fontsize',15)




x = LR(:,1)
y = LR(:,4)
err = LR(:,5)
figure(2)
errorbar(x,y,err,'color',[0 0 0],'linewidth',2)
hold on

x =SVR_RBF(:,1)
y = SVR_RBF(:,4)
err = SVR_RBF(:,5)
figure(2)
errorbar(x,y,err,'color',[0 0 1],'linewidth',2)

x =SVR_LR(:,1)
y = SVR_LR(:,4)
err = SVR_LR(:,5)
figure(2)
errorbar(x,y,err,'color',[1 0 0],'linewidth',2)
xlabel('Time [Sec]')
ylabel('RMSE [m]')
set(gca,'fontsize',15)


figure(3)
plot(x,LR(:,2),'k','linewidth',lw)
hold on
plot(x,SVR_RBF(:,2),'r','linewidth',lw)
plot(x,SVR_LR(:,2),'b','linewidth',lw)
h1 = plot(x,LR(:,2),'ko','linewidth',lw,'markersize',mk)
h2 = plot(x,SVR_RBF(:,2),'rx','linewidth',lw,'markersize',mk)
h3 = plot(x,SVR_LR(:,2),'bd','linewidth',lw,'markersize',mk)
xlabel('Time [Sec]')
ylabel('RMSE [m]')
title('k-Fold cross validation')
set(gca,'fontsize',15)
legend([h1, h2, h3],'LR', 'SVR-RBF', 'SVR-L')



figure(4)
plot(x,LR(:,4),'k','linewidth',lw)
hold on
plot(x,SVR_RBF(:,4),'r','linewidth',lw)
plot(x,SVR_LR(:,4),'b','linewidth',lw)
h1 = plot(x,LR(:,4),'ko','linewidth',lw,'markersize',mk)
h2 = plot(x,SVR_RBF(:,4),'rx','linewidth',lw,'markersize',mk)
h3 = plot(x,SVR_LR(:,4),'bd','linewidth',lw,'markersize',mk)
xlabel('Time [Sec]')
ylabel('RMSE [m]')
title('Leave-one-out cross validation')
set(gca,'fontsize',15)
legend([h1, h2, h3],'LR', 'SVR-RBF', 'SVR-L')



