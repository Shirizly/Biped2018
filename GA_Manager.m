addpath(genpath('MOGA_runs_scripts_and_functions'),genpath('Aux functions in use'),genpath('Stochastic terrains'),genpath('Results analysis'),...
    genpath('results 0921'));
% 
% clear all
% clc

%%
if exist('GA','var')
    %%
GA = GA.Run();
GA.Plot('Fit');
end

%% Comparison parameters:
% Controller Structure Selection
whichCPGp{1} = '6N_tagaLike_2Ank_torques_symm_feedback_eq'; 
whichCPGp{2} = 'ConSpitz';
whichCPGp{3} = 'ConSpitz_eq';
whichCPGp{4} = '2_level_CPG';
whichCPGp{5} = '2_level_CPG_eq';

CPGlist = [3];
adaptation = 0;

% Number of GA runs:
nRuns = 10;

% Evulotionary parameters:
MoogaGen = 20;
MoogaPop = 1000;

% Stochastic terrain parameters:
% setup for StoFit:
TerVarS = -6.5;
TerVarE = -4.5;
nseg = 5;
xend = 20;
nTerForSto = 10;
TerFileName = 'sto_ter_array_0508_';

%% Background setup

trainingDataFile = 'MatsRandomRes_6N_general_TrainingSet.mat'; % training data file name, for rea's reseacrch, hybrid GA/NN
whichGA_Case = '_GA_only_Feedback';


nStruct = length(CPGlist);
for i=1:nStruct
    whichCPGv{i} = whichCPGp{CPGlist(i)};
end

FileName_date = ['11_15_'];

failure = [];



for j = 1:nStruct
    whichCPG = whichCPGv{j}
    if adaptation == 1
        whichCPG = [whichCPG '_adaptation'];
        FileIn = [whichCPGv{j} '_09_21_' num2str(i) '.mat'];
    else
        FileIn = [];
    end
    generate_GenomeFile(whichCPG); % this functions holds definitions for the genetic sequenceFileName_prefix = 'AS_';
    FileName_start = [whichCPG,'_'];
    try
    for i = 1:nRuns
%     for i=1:nRuns  % running every method several times for stastistical analysis
        GA = MOOGA(MoogaGen,MoogaPop); % (generations,population size)
        FileName_extra1 = [num2str(i)];
        FileName_extra2 = [];
        GA.FileOut = [FileName_start,FileName_date,...
            FileName_extra1,FileName_extra2,'.mat'];
        GA.FileIn = FileIn;
        
        % Stochastic terrain setup
        GA.TerVarS = TerVarS;
        GA.TerVarE = TerVarE;
        GA.nseg = nseg;
        GA.xend = xend;
        GA.nTerForSto = nTerForSto;
        GA.TerFileName = TerFileName;
%         if j == 1
%             GA = GA.SetFittest(15,60,5);
%         else
%             GA = GA.SetFittest(15,15,5);
%         end
        
        %run
        MOGA_Run(GA,whichCPG,whichGA_Case,i)
    end
    catch err
        failure{end+1} = err;
        disp(err);
        save('failureslog.mat','failure');
    end
end