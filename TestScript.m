

%%
clear all; close all; clc;


% % % % % Tagalike CPG case study:
% whichCPG = '6N_TagaLike_symm'; % CPG case study name
% trainingDataFile = 'MatsRandomRes_6N_TagaLike_TrainingSet_2.mat'; % training data file name
% generate_GenomeFile('6N_tagaLike_2Ank_torques_symm_feedback'); %(enable CPG feedback)
% % generate_GenomeFile('6N_tagaLike_2Ank_torques_symm');          %(disable CPG feedback)

% % % % % general Matsuoka CPG case study:
whichCPG = '6N_general'; % CPG case study name
trainingDataFile = 'MatsRandomRes_6N_general_TrainingSet.mat'; % training data file name, for rea's reseacrch, hybrid GA/NN
generate_GenomeFile('6N_general_2Ank_torques'); % this functions holds definitions for the genetic sequence

%% GA only (No Feedback)
whichGA_Case = '_GA_only';

%% GA + NN (No Feedback)
whichGA_Case = '_NN_classi_only';

%% GA + Rescale (No Feedback)
whichGA_Case = '_rescale_only';

%% GA+classifier+rescale (No Feedback)
generate_GenomeFile('6N_tagaLike_2Ank_torques_symm_with_rescale')
whichGA_Case = '_NNclassi_and_rescale';

%% GA only (With Feedback)
whichGA_Case = '_GA_only_Feedback';

%%  GA + NN (With Feedback)
whichGA_Case = '_NN_classi_only_Feedback';

%% Begin Simulations:
for i=1:2  % running every method several times for stastistical analysis
    MOGA_Matsuoka_Run(whichCPG,whichGA_Case,trainingDataFile,[])
end
