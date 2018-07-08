function dx=matsuoka_der(t,x,a,b,c,tau,T)
% syms a b c tau T

A = [-1/tau 0 -b/tau 0;
    0 -1/tau 0 -b/tau;
    0 0 -1/T 0;
    0 0 0 -1/T];
B = [0 -a/tau;
    -a/tau 0;
    1/T 0;
    0 1/T];
y = max(x(1:2),0);
C = [ones(2,1);zeros(2,1)];
dx = A*x+B*y+c/tau*C;
end