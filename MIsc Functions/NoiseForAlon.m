clc;clear all;%close all

dt = 1/1000;
tf = 10;


t = 0:dt:tf;

f0 = 0.1;
f1 = 2;


y =t.*randn(1,length(t));
subplot 211
% plot(t,y)
hold on

s0 = 20;
sf = 1;
N = s0;

smooth_len = length(t)-N;

for i = (N+1):length(t)-(N+1)
    
smooth_len(i-N) = floor(s0 - s0*i/length(t) + sf );

t_tmp = t(i-N:i+N);
y_tmp = smooth( y(i-N:i+N) ,smooth_len(i-N));
y1(i-N) = y_tmp(N);
t1(i-N) = t(i);


end


plot(t1(1:2000),y1(1:2000),'r')
% title('Increasingly stochastic terrain function')
% xlabel('distance [m]')
% ylabel('height')
subplot 212
plot(t1(1:2000),smooth_len(1:2000))
% title('LPF frequency')
% xlabel('distance [m]')
% ylabel('frequency [Hz]')
hold on