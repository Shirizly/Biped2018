tv = 0:0.01:10;
w = 2*pi;
sig = 0.1;

syms t

R = [cos(w*t) sin(w*t);-sin(w*t) cos(w*t)];
s = exp(sig*t);

x0 = [1;0];
x = s*R*x0;
x=subs(x,t,tv);

plot(x(1,:),x(2,:))
axis equal
grid on