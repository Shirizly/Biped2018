clear all
load('LC_analysis')

for i=1:size(eigcell,2)
    eigcell{i}
    lpmax(i) = max(abs(real(eigcell{i})));
end


figure
plot(av1,lpmax,'lineWidth',2)
axis([av1(1),av1(end),0,1]);
ylabel('max | \lambda |');
xlabel('a - mutual inhibition parameter');
title('Max Eigen-value of the linearized Poincare map vs. the Mutual Inhibition Parameter')


figure
plot(av1,freqv)
ylabel('freq [Hz]');
xlabel('a - Mutual Inhibition Parameter');
title('Frequency of the LC period vs. the Mutual inhibition Parameter')