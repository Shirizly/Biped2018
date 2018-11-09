addpath(genpath('Results analysis'));
addpath(genpath('results 0921'));

%% Running the resulting controller of the GA process through a simulation

% first load GA reult .mat
seqs = GA.Seqs;
fits = GA.Fit;

GAend = GA.Progress;

% switch size(GA.Seqs,2)
% %     case 10
% %         start_pulses = 8; % for CPG
% %         np = 2;  
% %     case 19
% %         start_pulses = 8; % for CPG2
% %         np = 4;
% %         
% %     case 20
% %         start_pulses = 9; % for CPG3
% %         np = 4;
% end

% Taking the top controller for the i-st FF:
ind = find(fits(:,6,GAend)==max(fits(:,6,GAend)));
seq = seqs(ind,:,GAend)
fit = fits(ind,:,GAend)

%%
amp = seq(start_pulses:start_pulses+np-1).';

start = seq(start_pulses+np:start_pulses+2*np-1).';

stop = seq(start_pulses+np:start_pulses+2*np-1).'+seq(start_pulses+2*np:start_pulses+3*np-1).';

pulses = [amp,start,stop]
%% Taking the top controller for the second FF:
seq = seqs(2,:,GAend)
fit = fits(2,:,GAend)

amp = seq(start_pulses:start_pulses+np-1).';

start = seq(start_pulses+np:start_pulses+2*np-1).';

stop = seq(start_pulses+np:start_pulses+2*np-1).'+seq(start_pulses+2*np:start_pulses+3*np-1).';

pulses = [amp,start,stop]

%% Taking a balance controller:
nm = 3;
seq = seqs(nm,:,GAend)
fit = fits(nm,:,GAend)

amp = seq(start_pulses:start_pulses+np-1).';

start = seq(start_pulses+np:start_pulses+2*np-1).';

stop = seq(start_pulses+np:start_pulses+2*np-1).'+seq(start_pulses+2*np:start_pulses+3*np-1).';

pulses = [amp,start,stop]

%% Graphing the pulses:

phi = 0:0.0001:1;
p = zeros(np,size(phi,2));
on = ones(size(p));

p(phi>=start) = on(phi>=start);
p(phi>stop) = p(phi>stop)-on(phi>stop);
pul = [[amp(1:np/2).' zeros(1,np/2)] * p;[zeros(1,np/2) amp(np/2+1:np).'] * p];

figure
hold on
plot(phi,pul)
title('torque pulses of the CPG controller')
xlabel('Phase')
ylabel('magnitude [NM]')
legend('Ankle','Hip')
grid on
hold off


%% The simulation

ter_type = 1;

% sim = tryToWalk4(seq,ter_type);
sim = tryToWalkMatsuoka(seq,ter_type);

% ter = 5
% sim = tryToWalk2(seq,ter_type,ppv{ter},dppv{ter},whichCPG);


%% The simulation


sim = tryToWalk3(sim);



% %% low ankle torque:
% ind = find(seqs(:,4,end)<0.01);
% 
% ind2 = find(fits(ind,5
