
figure
hold on
plot(av,max(abs(lpv)),'lineWidth',2)
plot(av,max(abs(lpvf(:,1:998))),'lineWidth',2)
axis([floor(av(1)),floor(av(end))+1,0,5]);
line([floor(av(1)),floor(av(end))+1],[1,1])
ylabel('Max | \lambda |');
xlabel('Mutual Inhibition Parameter');
title('Max Eigen-value of the linearized Poincare map vs. the Mutual Inhibition Parameter')
hold off