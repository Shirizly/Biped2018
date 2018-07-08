tau = 0.0241*5;
beta = 18.9771;
tin = 4.7941;
a = 1.9531; %w in the controller

T = tau*12;

wn = 1/T*(sqrt(((tau+T)*beta-tau*a)/(tau*a)))
%%


MO = ControllerMatsuoka(tau,beta,tin,a);
tspan = 0:0.001:20;

[T,X] =ode45(@MO.Derivative,tspan,MO.IC);
y = X(:,1)-X(:,2);

% plot(T,X(:,1:2));
figure
plot(T,y);
grid on