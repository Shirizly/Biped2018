function [pp,dpp] = Sto_Ter_Gen2(tvar,xend)


% close all
% clear 
% clc

dt = 1/1000;
Fs = 1/dt;
tf = xend+dt;

t = 0:dt:tf;





y =2*(0.5*ones(1,length(t)) - rand(1,length(t)));

% plot(t,y)

% for i = 1:20 %moving average filter, discarded for butterworth filter
%     yf = (yf+[yf(2:end),yf(end)])/2;
% end
%%
[z,p,k] = butter(4,2*dt); %the cutoff is the number divided by 2
sos = zp2sos(z,p,k);
% fvtool(sos,'Analysis','freq','Fs',Fs) %for seeing the Bode of the filter
yf = sosfilt(sos,y);
%% match target variance
yvar = 1/size(yf,2)*sum((yf-mean(yf)*ones(1,length(yf))).^2);
yf = yf/sqrt(yvar/tvar);
var = 1/size(yf,2)*sum((yf-mean(yf)*ones(1,length(yf))).^2)
% figure
% plot(t,yf,'lineWidth',2)
% grid on
% xlabel('X [m]');
% ylabel('Y [m]');
% title(['Stochastically generated terrain (BL white noise), var = ' num2str(var)])
% axis([0 3 -0.2 0.2])
%% make the spline
n=10;
yf = [5 zeros(1,n) yf yf(end)*ones(1,n) 5];
t = [t(1)-(n+1)*100*dt:100*dt:t(1)-100*dt t t(end)+100*dt:100*dt:t(end)+(n+1)*100*dt];
pp = spline(t,yf);
dpp = mkpp(t,repmat(pp.order-1:-1:1,pp.dim*pp.pieces,1).*pp.coefs(:,1:pp.order-1),pp.dim);
dyf = ppval(dpp,t);
yr = dyf(n:end-n);
var = 1/size(yr,2)*sum((yr-mean(yr)*ones(1,length(yr))).^2);
% figure
% plot(t(n:end-n),180/pi*dyf(n:end-n))
% grid on
% xlabel('X [m]');
% ylabel('Y [m]');
% title('Stochastically generated terrain')

% axis equal
% Plotting the graphs
% figure
% plot(t,yf)
% axis equal
% grid on
% title(['stochastic terrain with Var =' num2str(var)]);
% xlabel('Distance [m]')
% ylabel('Height [m]')

% figure % This is a spectral analysis of the filtered/unfiltered terrain

L = length(y);
Y = fft(y);
Yf = fft(yf);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);

Pf2 = abs(Yf/L);
Pf1 = Pf2(1:L/2+1);

% P1 = P1*(Pf1(1)/P1(1));

f = Fs*(0:(L/2))/L;
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

end
