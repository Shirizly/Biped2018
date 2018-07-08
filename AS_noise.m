% plot(yout{1}.Values)
%%
clear
clc
amp = 10;
maxSl = 0.3;

dx = 1;
xend = 10;
N = xend/dx+1;
% X = (0:dt:tend);

% X = (randn(N,1)).';
% X = sort(X);
% X = X - min(X);
% X = X.*(xend/max(X));
X(1) = 0;
Y(1) = 0;
for i=2:N
    dx = min(sqrt(1/i)*abs(rand()),0.5);
    X(i) = X(i-1) + dx;
    dy = max(min(rand()-0.5,maxSl*dx),-maxSl*dx);
    Y(i) = Y(i-1) + dy;
end
X = X.*(xend/max(X));
dt = atan2(Y(end)+rand()-0.5,X(end));
R = [cos(dt) sin(dt);-sin(dt) cos(dt)];

for i=2:size(Y,2)
    XY = R*[X(i);Y(i)];
    X(i) = XY(1);
    Y(i) = XY(2);
end

% Y = randn(size(X,2),1).';
% Y = min(Y,1);
% Y = max(Y,-1);

% X = ds*yout{1}.Values.Time.';
% Y = yout{1}.Values.Data.';

% figure
% hold on
% plot(0:10,X)
% plot(X,Y)
%%

for j = 1:5
    X = [X,X(end)];
    Y = [Y,Y(end)];
    for i=1:(size(X,2)-1)
        XX(2*i) = (X(i)+ X(i+1))/2;
        XX(2*i-1) = X(i);
        YY(2*i-1) = Y(i);
        YY(2*i) = (Y(i)+ Y(i+1))/2;
    end
    X = XX(1:end-1);
    Y = YY(1:end-1);
    

    
end
% figure
% plot(X,Y)
%%
Lin = 1:1:(size(X,2));
T = Y.*sqrt(Lin)/(amp);
% figure
% plot(X,T)
X1 = 0:0.01:X(end);
T1 = spline(X,T,X1);



figure
hold on
plot(X1,T1)
title('Increasingly stochastic terrain function')
xlabel('distance [m]')
ylabel('height [m]')
axis equal
hold off
%%
% % Xn = ds*(0:dt*2:tend);
% Yn = smooth(Y,3).';
% Yn = spline(X,Yn,X1)/2;
% Tn = Yn.*(X1)/(amp);
% 
% XO = ones(size(X1,2),1).'*X(end);
% XL =XO-logspace(log(X(end)-1)/log(10),0,size(X1,2))+1;
% 
% figure
% hold on
% plot(XL,Tn)
% axis equal
% hold off