%% Printing fronts
clc

Fi = 1;
% gen = [30];


figure
hold on
grid on

% filename = {'AS_2_level_CPG2_withFA_02_20_03_18_30Gen_100pop__GA_only_Feedback',
cols = [1 0 0; 0 0 1];
contyp = {'2-lvl','6N-Taga'};
gen = {10,20};   

    uiopen('mat')
    
    Seqs = GA.Seqs;
    Fits = GA.Fit;
       
for j = 1:2
    % first load GA reult .mat
%     clear GA;

    
%   cols = rand(length(gen),3);
    

    
%     for i = 1:length(gen)
        fit = Fits(:,GA.FitIDs,gen{j});
        col = cols(j,:);
        GA.pareto_plot(GA,fit,Fi,col)
%         leg{i} = ['generation ' num2str(gen(i))];
%     end
%     legend(leg)
% contyp = input('controller type ','s');
% leg{j} = ['Con - ' contyp{j}];
leg{j} = ['gen - ' num2str(gen{j})];
end
legend(leg)
title('1st pareto front comparison - 30 gen,100 pop')
xlabel('Velocity FF');
ylabel('Energy Efficiency FF');
hold off



