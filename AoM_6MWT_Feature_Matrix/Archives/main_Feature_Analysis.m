clc
clear all
close all

tbl = import_Feature_Matrix('./MWT6_AoM_Feature_Matrix360.csv')

AoM = tbl{:,9:end}
corr(AoM)

% corr > 0.95
corr(AoM(:,2),AoM(:,10))
figure(1)
plot(AoM(:,2),AoM(:,10),'ko')

corr(AoM(:,6),AoM(:,9))
figure(2)
plot(AoM(:,6),AoM(:,9),'ko')

corr(AoM(:,6),AoM(:,11))
figure(3)
plot(AoM(:,6),AoM(:,11),'ko')

corr(AoM(:,6),AoM(:,12))
figure(4)
plot(AoM(:,6),AoM(:,12),'ko')














