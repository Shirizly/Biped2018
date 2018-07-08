clear all


%%
load('LC_analysis_for_comparison.mat')
tau = 0.0598;
T = 12*tau;
b = 9.8715;
c = 1.1603;

for run = 1:length(av)
    disp(run)
    xLC = xLCv{run};
    a = av(run);
    Period = 1/freq(run);
    IClim = xLC(1,:).';
    PMeps = 1e-3;
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
    
    
    t = 0:0.001:2*Period;
    EE=1e-9;
    options = odeset('Events',@ZC,'RelTol',EE,'AbsTol',EE);
    
    for d = Coords
        sIC = dIC(:,d);
        
        [t,x,te,ye,ie]=ode45(@matsuoka_der,t,sIC,options,a,b,c,tau,T);
%         figure
%         plot(t,x(:,1:2))
        LCend = find(ie(2:end)==5,1,'first')+1;
        dICp(:,d) = ye(LCend,:).';
    end
    
    % Calculate deviation
    DeltaM = dICp - IC;
    
    DP = 1/PMeps*(DeltaM);
    [EigVec,EigVal] = eig(DP,'nobalance');
    EigVal = diag(EigVal);
%     disp(EigVal)
%     lpv(:,run)
    lnv(:,run) = EigVal;
end

%%
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