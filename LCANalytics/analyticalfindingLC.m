syms t x tau T a b
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

A = AB{3}
eat = simplify(expm(A*t))
% 
% out = int(expm(A*(t-x)),x)
% input = subs(out,x,t)-subs(out,x,0)
in = simplify(inv(A)*expm(A*t))