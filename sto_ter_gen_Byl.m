% function [pp,dpp] = sto_ter_gen(n,slo_lim,xend)
% this function receives: n- number of changes in slope over the terrain length, slo_lim - limit on
% absolute value of slope, xend - length of terrain's X distance
% the function returns a random stochastic terrain (and its derivative) as
% piecewise polynomial structures.
% if you want constant stochastics, rather than increasing, set rr in both
% for loops to 1 (this is the parameter that deals with the growth is
% randomness)

clear all

n = 15;
tvar = 1E-4;
xend = 3;


xbase = 0.15; % this parameter defines the relative influence of the randomness in the horizontal lengths of the segments
% basically, each segment has a length which is the sum of xbase and a
% random number (between 0 and rr, which is linearly decreasing from 1).
% the actual segment length is different, since the x vector is later
% normalized to reach xend.


% x vector
xp(1) = 0;

for i = 2:n % creating the x vector, with stoch. varying horizontal distance
    
    xp(i) = xp(i-1)+(xbase+0.05*rand(1));   % each segment has a base length and a random addition
end

xp = xp*xend/xp(end);  % stretching the x vector values to end at xend


% y vector
yp(1) = 0;
avg_slo = 0;
slo = 0;
sig = 1.5*pi/180;
lim = 6*sig;
for i = 2:n % creating the y vector, with stoch. varying slope, limited by slo_lim, with a small bias to negate past slope
    dx = (xp(i)-xp(i-1));  
    rr = max(min(sig*randn(1),lim),-lim) % random number between -1 and 1
   slo(i) = 180/pi*rr;
    yp(i) = rr*dx+yp(i-1);
end
%match target variance
yvar = (1/n)*sum((yp-mean(yp)*ones(1,length(yp))).^2);
% yf = yp/sqrt(yvar/tvar);
% yvar = (1/n)*sum((yf-mean(yf)*ones(1,length(yf))).^2)
figure
plot(xp,yp,'-','lineWidth',2)
grid on
axis([0 3 -0.2 0.2])
xlabel('X [m]');
ylabel('Y [m]');
title(['Stochastically generated terrain (BYL), var= ' num2str(yvar)])
figure


plot(xp(2:end),slo(2:end),'-','lineWidth',2)
grid on
%%
% %temporary solution for good controllers, and controllers walking backwards
% xp = [-50, xp];
% xp(end+1) = xp(end)+10;
% 
% yp = [0,yp];
% yp(end+1) = yp(end)+10; % 0.5 slope for the extra part, for stopping good controllers from continuing beyind the end of the terrain
% 
% 
% % Adding midpoints for spline purposes
% X = xp;
% Y = yp;
% for j = 1:3 % Adding midpoints to the curve to make the spline more similar to it
%     for i=1:(size(X,2)-1) %For each point, it itself and the midpoint between it and the next are added to the new vector
%         XX(2*i) = (X(i)+ X(i+1))/2;
%         XX(2*i-1) = X(i);
%         YY(2*i-1) = Y(i);
%         YY(2*i) = (Y(i)+ Y(i+1))/2;
%     end
% 
% X = [XX X(end)]; % In the end of each run over the vectors, the end point needs to be added as well
% Y = [YY Y(end)];
% end
% 
% 
% 
% 
% pp = spline(X,Y);
% dpp = mkpp(X,repmat(pp.order-1:-1:1,pp.dim*pp.pieces,1).*pp.coefs(:,1:pp.order-1),pp.dim);
% 
% % Plot of the terrain:
% figure
% x = 0:0.01:xend;
% % x = min(X):0.01:max(X);
% y = spline(X,Y,x);
% plot(x,y,'LineWidth',2)
% xlabel('X [m]');
% ylabel('Y [m]');
% title('Stochastically generated terrain')
% axis equal
% grid on
% end