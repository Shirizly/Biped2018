% function [pp,dpp] = sto_ter_gen(n,slo_lim,xend)
% this function receives: n- number of changes in slope over the terrain length, slo_lim - limit on
% absolute value of slope, xend - length of terrain's X distance
% the function returns a random stochastic terrain (and its derivative) as
% piecewise polynomial structures.
% if you want constant stochastics, rather than increasing, set rr in both
% for loops to 1 (this is the parameter that deals with the growth is
% randomness)

clear all
xend = 10;
n = 10;
tvar = 1;
slo_lim = 0.3; %defines the absolute value of the largest slope allowed


xbase = 0.5; % this parameter defines the relative influence of the randomness in the horizontal lengths of the segments
% basically, each segment has a length which is the sum of xbase and a
% random number (between 0 and rr, which is linearly decreasing from 1).
% the actual segment length is different, since the x vector is later
% normalized to reach xend.


% x vector
xp(1) = 0;

for i = 2:n % creating the x vector, with stoch. varying horizontal distance
    
    xp(i) = xp(i-1)+(xbase+rand(1));   % each segment has a base length and a random addition
end

% xc = xc-ones(1,length(xc))*xc(1); % adjusting the x vector values to start at zero
xp = xp*xend/xp(end); % stretching the x vector values to end at xend


% y vector
yp(1) = 0;
avg_slo = 0;
slo = 0;
for i = 2:n % creating the y vector, with stoch. varying slope, limited by slo_lim, with a small bias to negate past slope
    rr = (1/(n))*(i-2); % the gradual growth function
    
    dx = (xp(i)-xp(i-1));  
    rr = (1-2*rand(1)); % random number between -1 and 1
   
    yp(i) = slo_lim*rr*dx+yp(i-1);
end
%match target variance
x = 0:0.01:xend;
y = spline(xp,yp,x);
yvar = sum((y-mean(y)*ones(1,length(y))).^2);
yf = y/sqrt(yvar/tvar);
var = sum((yf-mean(yf)*ones(1,length(yf))).^2)
plot(x,yf)
xlabel('X [m]');
ylabel('Y [m]');
title('Stochastically generated terrain')
axis equal
%%
%temporary solution for good controllers, and controllers walking backwards
xp = [-50, xp];
xp(end+1) = xp(end)+10;

yp = [0,yp];
yp(end+1) = yp(end)+10; % 0.5 slope for the extra part, for stopping good controllers from continuing beyind the end of the terrain


% Adding midpoints for spline purposes
X = xp;
Y = yp;
for j = 1:3 % Adding midpoints to the curve to make the spline more similar to it
    for i=1:(size(X,2)-1) %For each point, it itself and the midpoint between it and the next are added to the new vector
        XX(2*i) = (X(i)+ X(i+1))/2;
        XX(2*i-1) = X(i);
        YY(2*i-1) = Y(i);
        YY(2*i) = (Y(i)+ Y(i+1))/2;
    end

X = [XX X(end)]; % In the end of each run over the vectors, the end point needs to be added as well
Y = [YY Y(end)];
end




pp = spline(X,Y);
dpp = mkpp(X,repmat(pp.order-1:-1:1,pp.dim*pp.pieces,1).*pp.coefs(:,1:pp.order-1),pp.dim);

% Plot of the terrain:
figure
x = 0:0.01:xend;
% x = min(X):0.01:max(X);
y = spline(X,Y,x);
plot(x,y,'LineWidth',2)
xlabel('X [m]');
ylabel('Y [m]');
title('Stochastically generated terrain')
axis equal
grid on
end


