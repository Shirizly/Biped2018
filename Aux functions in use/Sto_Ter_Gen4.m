function [pp,dpp] = Sto_Ter_Gen4(tvarstart,tvarend,nseg,xend)

%% part for debugging
% close all
% clear 
% clc
% 
% tvarstart = -5;
% tvarend = -4;
% nseg = 5;
% xend = 25;
flag = 1;
%% creating the increasing variance white noise

varv = logspace(tvarstart,tvarend,nseg)

dt = 1/1000;
Fs = 1/dt;
tend = xend+dt;
varnormsr = sqrt(524*12);
t = 0:dt:tend;
while flag

y = varnormsr * (0.5*ones(1,length(t)) - rand(1,length(t)));
% var = 1/size(y,2)*sum((y-mean(y)*ones(1,length(y))).^2);

segstr = floor(linspace(1,length(y),nseg+1));

for i = 1:nseg
    y(segstr(i):segstr(i+1)) = y(segstr(i):segstr(i+1))*sqrt(varv(i));
    var(i) = 1/524*1/size(y(segstr(i):segstr(i+1)),2)*sum((y(segstr(i):segstr(i+1))-mean(y(segstr(i):segstr(i+1)))*ones(1,length(y(segstr(i):segstr(i+1))))).^2);
end

% figure
% plot(t,y)


[z,p,k] = butter(4,2*dt); %the cutoff is the number divided by 2
sos = zp2sos(z,p,k);
% fvtool(sos,'Analysis','freq','Fs',Fs) %for seeing the Bode of the filter
yf = sosfilt(sos,y);

%% match target variance
yvar = 1/size(yf,2)*sum((yf-mean(yf)*ones(1,length(yf))).^2);

for i = 1:nseg
%      var1(i) = 1/size(yf(segstr(i):segstr(i+1)),2)*sum((yf(segstr(i):segstr(i+1))-mean(yf(segstr(i):segstr(i+1)))*ones(1,length(yf(segstr(i):segstr(i+1))))).^2);
%      yf(segstr(i):segstr(i+1)) = yf(segstr(i):segstr(i+1))*sqrt(varv(i)/var1(i));
     var2(i) = 1/size(yf(segstr(i):segstr(i+1)),2)*sum((yf(segstr(i):segstr(i+1))-mean(yf(segstr(i):segstr(i+1)))*ones(1,length(yf(segstr(i):segstr(i+1))))).^2);

     diffv(i) = (var2(i)-var(i))/var(i);
end
diffv
flag = max(abs(diffv))>0.2;
    
% for i = 1:nseg % Fitting the variance after the filter, creates jumps
%     yf(segstr(i):segstr(i+1)) = yf(segstr(i):segstr(i+1))*sqrt(varv(i));
%     var2(i) = 1/size(yf(segstr(i):segstr(i+1)),2)*sum((yf(segstr(i):segstr(i+1))-mean(yf(segstr(i):segstr(i+1)))*ones(1,length(yf(segstr(i):segstr(i+1))))).^2);
% end

% figure
% plot(t,yf,'lineWidth',2)
% grid on
% xlabel('X [m]');
% ylabel('Y [m]');
% title(['Stochastically generated terrain (BL white noise), var range = ' num2str(var2(1)) ' , ' num2str(var2(nseg))])
% % axis([0 xend -0.2 0.2])
% axis equal

%% make the spline
n=10;
yf = [5 zeros(1,n) yf yf(end)*ones(1,n) 5];
t = [t(1)-(n+1)*100*dt:100*dt:t(1)-100*dt t t(end)+100*dt:100*dt:t(end)+(n+1)*100*dt];
pp = spline(t,yf);
dpp = mkpp(t,repmat(pp.order-1:-1:1,pp.dim*pp.pieces,1).*pp.coefs(:,1:pp.order-1),pp.dim);


%% Checking the slopes
t = 0:0.001:xend;
dyf = 180/pi*ppval(dpp,t);
max(abs(dyf))
flag = max(abs(dyf))>10 || flag;
end
figure
t = 0:0.001:xend;
dyf = 180/pi*ppval(dpp,t);
plot(t,dyf)
% yr = dyf(n:end-n);
% var = 1/size(yr,2)*sum((yr-mean(yr)*ones(1,length(yr))).^2);
% figure
% plot(t,180/pi*dyf)
% grid on
% xlabel('X [m]');
% ylabel('Y [m]');
% title('Slope')

% axis equal
% Plotting the graphs
% figure
% plot(t,yf)
% axis equal
% grid on
% title(['stochastic terrain with Var =' num2str(var)]);
% xlabel('Distance [m]')
% ylabel('Height [m]')

%% This is a spectral analysis of the filtered/unfiltered terrain

% figure 
% L = length(y);
% Y = fft(y);
% Yf = fft(yf);
% 
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% 
% Pf2 = abs(Yf/L);
% Pf1 = Pf2(1:L/2+1);

% P1 = P1*(Pf1(1)/P1(1));

% f = Fs*(0:(L/2))/L;
% figure
% hold on
% ne = floor(length(f)/20);
% f=f(1:ne);
% P1 = P1(1:ne);
% Pf1 = Pf1(1:ne);
% plot(f,P1)
% plot(f,Pf1)
% title('Amplitude spectrum of the terrain');
% xlabel('frequency [1/m]');
% ylabel('amplitude')
% legend('unfiltered','filtered')

% end
