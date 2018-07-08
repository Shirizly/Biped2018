clear all
%%
% load('LC_analysis_new.mat')

for i = 1:size(av,2)
    maxl = max(abs(lpv(:,i)));
    k1 = find(abs(lpv(:,i))>0.0001*maxl);
    lpnztemp = lpv(k1,i);
    lpnz(1:size(k1,1),i) = sort(lpnztemp,'ComparisonMethod','real');
       
end

%%
n=1:3 ;
start = 1;
% figure
% hold on
plot(real(lpnz(n,start:end)).',imag(lpnz(n,start:end)).','.','markerSize',10)
% xlabel('Real')
% ylabel('imaginary')
% title(['EigenValues of the Linearized poincare map, a>= ' num2str(av(start))])
% grid on
% axis equal
% tet=0:0.001:2*pi;
% plot(cos(tet),sin(tet))

% figure
% plot(av,real(lpnz),'.')
% figure
% plot(av,imag(lpnz),'.')

%%
