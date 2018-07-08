clear all
clc

%% find limit cycle;
% try 
tau = 0.0598;
T = 12*tau;
b = 9.8715;
c = 1.1603;
% a = 7.2618;
n=100;
av = linspace(1.2,6,n);
% limits (1.1,10.87)

x0 = [1 0 0 0];



dt = 0.01;
tend = 200;


EE=1e-10;
options = odeset('Events',@ZC,'RelTol',EE,'AbsTol',EE);


tol1 = 1e-8;

skip = 0;
avnil = [];
%% Defining the constant matrices/vectors

HB{1} = [1 0 0 0];
HB{2} = [0 1 0 0];
HB{3} = [1 -1 0 0];


B = c*[1/tau;1/tau;0;0];

filter = zeros(n,1);
%%
run = 0;
while run<length(av)
    run = run+1
    a = av(run);
    flag = 1;    
    tv = 0:dt:tend;
    tend1 = tend;
    x1 = x0;
    if run>1 % take initial conditions from the LC of the previous iteration, assuming close a values, should be helpful to convergence
        x1 = xLCv{run-1}(1,:);
    end
    
    % Running the simulation until convergence
    
    [t,x,te,ye,ie]=ode45(@matsuoka_der,tv,x1,options,a,b,c,tau,T);
    LCend = find(ie==5,1,'last');
    LCstartv = find(ie(1:LCend-1)==5);
    
    i=0;
    while flag
        LCstart = LCstartv(end-i);
        
        if max(abs(ye(LCstart,:)-ye(LCend,:)))<tol1
            flag = 0;
            yLCstart = ye(LCend,:);
            disp(i)
%             convest = max(abs(ye(LCstart,:)-ye(LCend,:)))
        else
            i = i+1;
        end
    end
        
        
%         if max(abs(ye(LCstart,:)-ye(LCend,:)))<tol1
%             flag = 0;
%         else
%             if tend1<2^10*tend
%                 convest = max(abs(ye(LCstart,:)-ye(LCend,:)))
%                 tv = tend1:dt:4*tend1;
%                 tend1 = tend1*4
%                 x1 = x(end,:);
%             else
%                 flag = 0;
%                 disp('converging takes longer than reasonable time!');
%                 skip = 1;
%             end
%         end
%     end
    
    nZC = LCend-LCstart-1;
%     % Dealing with skips in case of nonsufficient convergence of the LC
%     if skip == 1                
%         av = [av(1:run-1),av(run+1:end)];
%         run = run-1;
%         skip = 0;
%         continue;
%     end
    
    %% Building vectors defining the ZC

    for i = 1:nZC+2 
        ind(i) = find(abs(t-te(LCend-nZC+i-2))<dt/2);
        dir(i) = ie(LCend-nZC+i-2);
        xe(i,:) = ye(LCend-nZC+i-2,:);
    end
    
    %% removing 2 close ZC that are a result of grazing
%     dirstr = num2str(dir);
%     pat = {'2  1  2  1','1  2  1  2','3  4  3  4','4  3  4  3'};
%     for i = 1:4   
%     if contains(dirstr,pat{i})
%         in = strfind(dirstr,pat{i});
%         remove = [floor(in/3)+2,floor(in/3)+3];
%         dirstr(in+3:in+8) = [];
%         ind(remove) = [];
%         dir(remove) = [];
%         xe(remove,:) = [];
%         nZC = nZC - 2;
%         filter(run) = filter(run)+1;
%     end
%     end
    
    
    xLC = x(ind(1):ind(2+nZC),:);
    tLC = t(ind(1):ind(2+nZC))-ones(ind(2+nZC)-ind(1)+1,1)*t(ind(1));
    
    
%     ind2 = find(abs(t-te(end-nZC*2))<dt/2);     
%     xLC2 = x(ind2:ind(1),:);
%     tLC2 = t(ind2:ind(1))-ones(ind(1)-ind2+1,1)*t(ind2);
    
    ind = round(ind-ones(1,length(ind))*(ind(1)-1),0);


%%

period = tLC(end);
freq(run) = 1/period;
ind0 = ind(1);
ind = ind(2:nZC+2);

AB{1} = [-1/tau -a/tau -b/tau 0; %% both are positive
    -a/tau -1/tau 0 -b/tau;
    1/T 0 -1/T 0;
    0 1/T 0 -1/T];

AB{2} = [-1/tau 0 -b/tau 0; %% u2 is negative
    -a/tau -1/tau 0 -b/tau;
    1/T 0 -1/T 0;
    0 0 0 -1/T];

AB{4} = [-1/tau -a/tau -b/tau 0; %% u1 is negative
    0 -1/tau 0 -b/tau;
    0 0 -1/T 0;
    0 1/T 0 -1/T];

AB{3} = [-1/tau 0 -b/tau 0; %% both are negative
    0 -1/tau 0 -b/tau;
    0 0 -1/T 0;
    0 0 0 -1/T];

for i = 1:nZC+1
    if xLC(ind(i)-1,1)>0
        if xLC(ind(i)-1,2)>0
%             disp('both positive')
            A{i} = AB{1};
        else
%             disp('1 positive')
            A{i} = AB{2};
        end
    else if xLC(ind(i)-1,2)>0
%             disp('2 positive')
            A{i} = AB{4};
        else
%             disp('both negative')
            A{i} = AB{3};
        end
    end
    H{i} = HB{floor((dir(i+1)+1)/2)};
end

ts(1) = tLC(ind(1));
for i = 2:nZC+1
ts(i) = tLC(ind(i))-tLC(ind(i-1));
end

%%
P=1;
tp = [];
% xp(1,:) = xe(1,:);
for j = 1:nZC+1
    Jphi = expm(A{j}*ts(j));
    F = A{j}*xe(j+1,:).'+B;
    SM = eye(4)-(F*H{j})/(H{j}*F);
    P = SM*Jphi*P;
    
    if abs(H{j}*F) < 1e-2 % Checking for grazing
        avnil = [avnil a];
        tp = [tp te(j+1)];
        indgraz = j;
        disp(j)
    end
    
%     intJphi = integral(@(h)expm(A{j}.*h),0,ts(j),'ArrayValued',1);
%     xp(j+1,:) = (Jphi*xp(j,:).'+intJphi*B).';
end

lp = eig(P)

lpv(:,run)=lp;

%%
Period = 1/freq(run);
IClim = yLCstart.';
PMeps = 1e-6;
Coords = 1:4;
Ncoord = 4;

% Limit cycle initial conditions
IC = repmat(IClim, 1, Ncoord);

% Disturbed initial conditions
dIC = IC;
dICp = IC;
for d = Coords
    dIC(d,d) = dIC(d,d) + PMeps;
end


% Run the simulations


tn = 0:0.001:1.05*Period;
EE=1e-10;
options = odeset('Events',@ZC,'RelTol',EE,'AbsTol',EE);

for d = Coords
    sIC = dIC(:,d);
    
    [tn,x,te,ye,ie]=ode45(@matsuoka_der,tn,sIC,options,a,b,c,tau,T);
%             figure
%             plot(tn,x(:,1:2))
    LCend = find(ie==5,1,'last');
    dICp(:,d) = ye(LCend,:).';
end

% Calculate deviation
DeltaM = dICp - IC;

DP = 1/PMeps*(DeltaM);
[EigVec,EigVal] = eig(DP,'nobalance');
EigVal = diag(EigVal)
%     disp(EigVal)
%     lpv(:,run)
lnv(:,run) = EigVal;
%% Loading specific LC from file
% run = 49;
% xLC = xLCv{run};
% tLC = tLCv{run};
% nZC = nZCv{run};
% ind = indv{run};
% a = av(run)

%% Graphing the LC
% figure
% hold on
% plot(tLC,xLC(:,1:2),'lineWidth',3);
% if ~isempty(tp)
% plot(tp,zeros(size(tp)),'o')
% end
% % plot(tLC2,xLC2(:,1:2),'--','lineWidth',3);
% for i = 1:nZC+1
%     line([tLC(ind(i)),tLC(ind(i))],[-3 1],'color','black','LineStyle','--','lineWidth',2)
% end
% title(['Membrane potential of a Matsuoka oscillator - over a single LC period - a = ' num2str(a)])
% xlabel('time [sec]')
% ylabel('membrane potential');
% legend('u_1','u_2')
% grid on
% hold off

% figure
% plot(xLC(:,1),xLC(:,2))
% for i = 1:nZC
%     
% line([xLC(ind(i)),xLC(ind(i))],[-3 1],'color','black','LineStyle','--','lineWidth',2)
% end

%%
xLCv{run}= xLC;
tLCv{run} = tLC;
nZCv{run} = nZC;
indv{run} = ind;
end
%%

% figure
% plot(av,max(abs(lpv)),'lineWidth',2)
% axis([floor(av(1)),floor(av(end))+1,0,5]);
% line([floor(av(1)),floor(av(end))+1],[1,1])
% ylabel('Max | \lambda |');
% xlabel('Mutual Inhibition Parameter');
% title('Max Eigen-value of the linearized Poincare map vs. the Mutual Inhibition Parameter')

figure
hold on
plot(av,max(abs(lpv)),'lineWidth',2)
plot(av,max(abs(lnv)),'--','lineWidth',2)
axis([floor(av(1)),floor(av(end))+1,0,0.4]);
% line([floor(av(1)),floor(av(end))+1],[1,1])
ylabel('Max | \lambda |');
xlabel('Mutual Inhibition Parameter');
title('Max Eigen-value of the linearized Poincare map vs. the Mutual Inhibition Parameter')
legend('Analytical','Numerical')
grid on
hold off


figure
plot(av,freq)
ylabel('Frequency [Hz]');
xlabel('Mutual Inhibition Parameter');
title('Frequency of the LC period vs. the Mutual Inhibition Parameter')
%%

% save('LC_analysis_zoom2_filtered','av','b','c','tau','T','lpv','freq','xLCv','tLCv','nZCv','indv','filter')
% save('LC_analysis_zoom2_unfiltered','av','b','c','tau','T','lpv','freq')
save('LC_analysis_for_comparison_zoom','av','b','c','tau','T','lpv','freq','xLCv','tLCv','nZCv','indv','filter','lnv')
%%
% catch
%     disp('there was en error');
% end