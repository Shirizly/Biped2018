%% Running the resulting controller of the GA process through a simulation

% first load GA reult .mat
seqs = GA.Seqs;
fits = GA.Fit;
ppv = GA.ppv;
dppv = GA.dppv;

GAend = GA.Generations-2;
seqlen = size(GA.Seqs,2);
switch seqlen
    case 13
        start_pulses = 8; % for CPG
        np = 2;

    
    case 19
        start_pulses = 8; % for CPG2
        np = 4;
        
    case 20
        start_pulses = 9; % for CPG3
        np = 4;
end
%% Taking the top controller for the first FF:
seq = seqs(1,:,GAend);
fit = fits(1,:,GAend)


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
nm = 29;
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
ter_type = 2;
sim = tryToWalk(seq,ter_type);

%% Set the simulation for StoFit
    % Set up the simulations
    pp=ppv{GAend};
    dpp = dppv{GAend};
    
% Dur = 180;
    sim = tryToWalk(seq,ter_type,pp,dpp,0);
