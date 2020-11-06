clc
clear all
close all


Hz = 1000;
dt = 1/Hz;
t = [1:1:1000]*dt;

A = 10
R = 50
w = 100
y = A*sin(w*2*pi*t)+R*rand(1,length(t))

plot(t,y)

sample_ent = sampen(y,1,0.2)