addpath(genpath('MOGA_runs_scripts_and_functions'),genpath('Aux functions in use'));

clear all
clc


%% Comparison parameters:

whichCPGp{1} = '6N_tagaLike_2Ank_torques_symm_feedback_eq'; 
whichCPGp{2} = 'ConSpitz';
whichCPGp{3} = 'ConSpitz_eq';
whichCPGp{4} = '2_level_CPG';
whichCPGp{5} = '2_level_CPG_eq'; 

CPGlist = [3,2];
nRuns = 2;

% Evulotionary parameters:
MoogaGen = 5;
MoogaPop = 50;

% Stochastic terrain parameters:
% setup for StoFit:
TerVarS = -5;
TerVarE = -3;
nseg = 5;
xend = 25;
nTerForSto = 10;
TerFileName = 'sto_ter_array_0807_';

%% Background setup

trainingDataFile = 'MatsRandomRes_6N_general_TrainingSet.mat'; % training data file name, for rea's reseacrch, hybrid GA/NN
whichGA_Case = '_GA_only_Feedback';


nStruct = length(CPGlist);
for i=1:nStruct
    whichCPGv{i} = whichCPGp{CPGlist(i)};
end

FileName_date = [datestr(now,'mm_dd') '_'];




for j = 1:nStruct
    whichCPG = whichCPGv{j}
    generate_GenomeFile(whichCPG); % this functions holds definitions for the genetic sequenceFileName_prefix = 'AS_';
    FileName_start = [whichCPG,'_'];
    
    for i=1:nRuns  % running every method several times for stastistical analysis
        GA = MOOGA(MoogaGen,MoogaPop); % (generations,population size)
        FileName_extra1 = [num2str(i)];
        FileName_extra2 = [];
        GA.FileOut = [FileName_start,FileName_date,...
            FileName_extra1,FileName_extra2,'.mat'];
        
        GA.TerVarS = TerVarS;
        GA.TerVarE = GA.TerVarE;
        GA.nseg = nseg;
        GA.xend = xend;
        GA.nTerForSto = nTerForSto;
        GA.TerFileName = TerFileName;
        
        MOGA_Run(GA,whichCPG,whichGA_Case,i)
    end
end