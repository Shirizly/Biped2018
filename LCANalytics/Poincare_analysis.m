% %% Stability Analysis
% syms a b c tau T
% X = sym('x',[4,1])
% 
% H1 = [1 0 0 0];
% H2 = [0 1 0 0];
% A1 = [-1/tau -a/tau -b/tau 0;
%     -a/tau -1/tau 0 -b/tau;
%     1/T 0 -1/T 0;
%     0 1/T 0 -1/T];
% 
% A2 = [-1/tau 0 -b/tau 0;
%     -a/tau -1/tau 0 -b/tau;
%     1/T 0 -1/T 0;
%     0 0 0 -1/T];
% A3 = A1;
% 
% A4 = [-1/tau -a/tau -b/tau 0;
%     0 -1/tau 0 -b/tau;
%     0 0 -1/T 0;
%     0 1/T 0 -1/T];
% 
% B = c*[1/tau;1/tau;0;0];
% F1 = A1*X+B;
% F2 = A2*X+B;
% F3 = F1;
% F4 = A4*X+B;
% 
% ksi1 = simplify(H2/(H2*F1))
% ksi2 = simplify(H2/(H2*F2))
% ksi3 = simplify(H1/(H1*F3))
% ksi4 = simplify(H1/(H1*F4))


%%
clear all
clc

load('LC_1.mat')
ind1 = ind;
ind = ind(2:5);

HB{1} = [1 0 0 0];
HB{2} = [0 1 0 0];

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


B = c*[1/tau;1/tau;0;0];

for i = 1:4
    if xLC(ind(i)-1,1)>0
        if xLC(ind(i)-1,2)>0
            A{i} = AB{1};
        else
            A{i} = AB{2};
        end
    else if xLC(ind(i)-1,2)>0
            A{i} = AB{3};
        else
            A{i} = AB{4};
        end
    end
    if abs(xLC(ind(i),1))<abs(xLC(ind(i),2))
        H{i} = HB{1};
    else
        H{i} = HB{2};
    end
end



t(1) = tLC(ind(1));
t(2) = tLC(ind(2))-t(1);
t(3) = tLC(ind(3))-t(2)-t(1);
t(4) = tLC(end)-t(3)-t(2)-t(1);

P=1;

for j = 1:4
    Jphi = expm(A{j}*t(j))
    F = A{j}*xLC(ind(j),:).'+B
    SM = eye(4)-(F*H{j})/(H{j}*F)
    P = SM*Jphi*P
    lp = eig(P)
end

