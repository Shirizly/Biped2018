clear all

%%
childstart = GA.Fittest(1)+GA.Fittest(2)+1;
childend = GA.Population-GA.Fittest(3);
childnum = childend-childstart+1;

stt = 2;

for i = stt:GA.Progress
    Count(i) = length(find(GA.Fit(:,1,i)<1));
    Countchild(i) = length(find(GA.Fit(childstart:childend,1,i)<1));
end
ratio = Count/GA.Population;
ratiochild = Countchild./childnum;
ratioadult = (Count-Countchild)/(GA.Population-childnum);
plot(stt:GA.Progress,[ratio(stt:end);ratioadult(stt:end);ratiochild(stt:end)])