syms T tau a b c 

B = c*[1/tau;1/tau;0;0];
AB{1} = [-1/tau -a/tau -b/tau 0; %% both are positive
    -a/tau -1/tau 0 -b/tau;
    1/T 0 -1/T 0;
    0 1/T 0 -1/T];

AB{2} = [-1/tau 0 -b/tau 0; %% u2 is negative
    -a/tau -1/tau 0 -b/tau;
    1/T 0 -1/T 0;
    0 0 0 -1/T];

AB{3} = [-1/tau -a/tau -b/tau 0; %% u1 is negative
    0 -1/tau 0 -b/tau;
    0 0 -1/T 0;
    0 1/T 0 -1/T];

AB{4} = [-1/tau 0 -b/tau 0; %% both are negative
    0 -1/tau 0 -b/tau;
    0 0 -1/T 0;
    0 0 0 -1/T];

for i =1:4
  A =  AB{i}
% det(eye*l-AB{i})
[v,l] = jordan(A);
disp('eigenvectors:')
pretty(v)
disp('eigenvalues:')
pretty(l)
end