%% Printing impulses
clear all
clc

% Fi = 1;
% gen = [30];

uiopen('mat')


cols = [1 0 0; 0 0 1];    

seqs = GA.Seqs;
Fits = GA.Fit(:,GA.FitIDs,end);
fronts = GA.get_fronts(GA,Fits);

%%
    
switch size(GA.Seqs,2)
    case 19
        start_pulses = 8; % for CPG2
        
    case 20
        start_pulses = 9; % for CPG3
end    

figure
hold on
    
for i = 1:length(fronts{1})
    
       

seq = seqs(fronts{1}(i),:,end);

amp = seq(start_pulses:start_pulses+3).';

start = seq(start_pulses+4:start_pulses+7).';

dur = seq(start_pulses+8:start_pulses+11).';

impulse_hip = sort(abs(amp(3:4)).*dur(3:4));

impulse_ratio_hip(i) = impulse_hip(1)/impulse_hip(2);

impulse_ankle = sort(abs(amp(1:2)).*dur(1:2));

impulse_ratio_ankle(i) = impulse_ankle(1)/impulse_ankle(2);

ih(i,:)=impulse_hip;
ia(i,:)=impulse_ankle;

% leg{j} = ['gen - ' num2str(gen{j})];
end

% plot(impulse_ratio_hip,impulse_ratio_ankle,'*')
plot(ih(:,1),ih(:,2),'*')

% legend(leg)
title('impulse magnitude in the 1st pareto front - hip joint')
xlabel('smaller impulse');
ylabel('stronger impulse');
grid on

hold off


figure
plot(ia(:,1),ia(:,2),'*')

title('impulse magnitude in the 1st pareto front - ankle joint')
xlabel('smaller impulse');
ylabel('stronger impulse');
grid on

hold off