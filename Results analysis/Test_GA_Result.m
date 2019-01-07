% addpath(genpath('Results analysis'));
addpath(genpath('results 1115'));
%%
% clear all

%% Running the resulting controller of the GA process through a simulation
% generate_GenomeFile('6N_tagaLike_2Ank_torques_symm_feedback_eq')
generate_GenomeFile('ConSpitz_eq_adaptation')

% first load GA reult .mat
seqs = GA.Seqs;
fits = GA.Fit;

GAend = GA.Progress

%%
GAend = 1;
np = 3;
start_pulses = 2;


%% Taking the top controller for the i-st FF:
% choose top fitness:
ff = 5;

ind = find(fits(:,ff,GAend)==max(fits(:,ff,GAend)));
seq = seqs(ind,:,GAend)
fit = fits(ind,:,GAend)
% k_amp = seq(12);

amp0 = seq(start_pulses:3:start_pulses+np*3-1).';

start = seq(start_pulses+1:3:start_pulses+np*3).';

stop = seq(start_pulses+2:3:start_pulses+np*3+1).';
%%
slv = [5 0 -5]
for i = 1:3
% slope = slv(i);
% 
% amp = amp0 + slope*k_amp*[-1;1;1];
% 
% for j = 1:3
%     amp(j) = min(max(amp(j),-50),50);
% end
% 
% pulses = [amp,start,stop]
% 
% %% Graphing the pulses:
% 
% phi = 0:0.0001:1;
% p = zeros(np,size(phi,2));
% on = ones(size(p));
% 
% p(phi>=start) = on(phi>=start);
% p(phi>stop) = p(phi>stop)-on(phi>stop);
% switch np
%     case 3
%         pul = [[amp(1).' zeros(1,2)] * p;[zeros(1,1) amp(2:3).'] * p];
%     case 4
%         pul = [[amp(1:np/2).' zeros(1,np/2)] * p;[zeros(1,np/2) amp(np/2+1:np).'] * p];
% end
% figure
% hold on
% plot(phi,pul,'lineWidth',4)
% title(['Torque pulses of the CPG controller, slope = ' num2str(slope)])
% xlabel('Phase')
% ylabel('Magnitude [NM]')
% legend('Ankle','Hip')
% grid on
% hold off
% 
end
%% The simulation

ter_type = 2;

% sim = tryToWalk4(seq,ter_type);
% sim = tryToWalkMatsuoka(seq,ter_type);
ter = 1
% filename = 'rec2level3.gif';
% sim = tryToWalk4(seq,ter_type,ppv{ter},dppv{ter},filename);
% sim = tryToWalkMatsuoka(seq,ter_type,ppv{ter},dppv{ter},filename);

ter_type = 1;
% filename = 'rec2levelslow.gif';
% sim = tryToWalk4(seq,ter_type);
% sim = tryToWalkMatsuoka(seq,ter_type,1,filename);
%% Stick figure plot:
filename = [];
Dur = 20;
timestep = 0.01;
geneNum = GAend;
GenID = ind;
ter_type = 3;

switch ter_type
    case 1
        CBstick_Figure_plot(GA,geneNum, GenID, Dur, timestep, filename, ter_type);
        
    case 2
        load('sto_ter_array_0508_11','ppv','dppv')
        ter = 1;
        CBstick_Figure_plot(GA,geneNum, GenID, Dur, timestep, filename, ter_type,ppv{ter},dppv{ter});
        
    case 3
        slope = 5;
        CBstick_Figure_plot(GA,geneNum, GenID, Dur, timestep, filename, ter_type,slope);
end
