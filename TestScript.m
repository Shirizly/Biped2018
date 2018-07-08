addpath(genpath('MOGA_runs_scripts_and_functions'),genpath('Aux functions in use'));

%%
clear all;
% close all; 
clc;


% % % % % Tagalike CPG case study:
% whichCPG = '6N_TagaLike_symm'; % CPG case study name
% trainingDataFile = 'MatsRandomRes_6N_TagaLike_TrainingSet_2.mat'; % training data file name
% generate_GenomeFile('6N_tagaLike_2Ank_torques_symm_feedback'); %(enable CPG feedback)
% % generate_GenomeFile('6N_tagaLike_2Ank_torques_symm');          %(disable CPG feedback)

% % % % % general Matsuoka CPG case study:
whichCPG = '6N_tagaLike_2Ank_torques_symm_feedback'; % CPG case study name
trainingDataFile = 'MatsRandomRes_6N_general_TrainingSet.mat'; % training data file name, for rea's reseacrch, hybrid GA/NN
generate_GenomeFile('6N_tagaLike_2Ank_torques_symm_feedback'); % this functions holds definitions for the genetic sequence

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
for i=1:1  % running every method several times for stastistical analysis
    MOGA_Matsuoka_Run(whichCPG,whichGA_Case,trainingDataFile,[],i)
end

%%
% part 2:
clear all; 
% close all; 
clc;


% % % % % general Matsuoka CPG case study:
whichCPG = 'ConSpitz'; % CPG case study name
trainingDataFile = 'MatsRandomRes_6N_general_TrainingSet.mat'; % training data file name, for rea's reseacrch, hybrid GA/NN
generate_GenomeFile(whichCPG); % this functions holds definitions for the genetic sequence
startdate = datestr(now,'mm_dd');

%  GA only (With Feedback)
whichGA_Case = '_GA_only_Feedback';

%  Begin Simulations:
for i=1:10  % running every method several times for stastistical analysis
    MOGA_Matsuoka_Run(whichCPG,whichGA_Case,trainingDataFile,[],startdate,i)
end

% % % % % general Matsuoka CPG case study:
whichCPG = 'ConSpitz_eq'; % CPG case study name
trainingDataFile = 'MatsRandomRes_6N_general_TrainingSet.mat'; % training data file name, for rea's reseacrch, hybrid GA/NN
generate_GenomeFile(whichCPG); % this functions holds definitions for the genetic sequence


%  GA only (With Feedback)
whichGA_Case = '_GA_only_Feedback';

%  Begin Simulations:
for i=1:10  % running every method several times for stastistical analysis
    MOGA_Matsuoka_Run(whichCPG,whichGA_Case,trainingDataFile,[],startdate,i)
end
%%
% clear all; 
% % close all; 
% clc;




% % % % % general Matsuoka CPG case study:
whichCPG = '2_level_CPG3'; % CPG case study name
trainingDataFile = 'MatsRandomRes_6N_general_TrainingSet.mat'; % training data file name, for rea's reseacrch, hybrid GA/NN
generate_GenomeFile('2_level_CPG3'); % this functions holds definitions for the genetic sequence

%  GA only (With Feedback)
whichGA_Case = '_GA_only_Feedback';

%  Begin Simulations:
for i=1:10  % running every method several times for stastistical analysis
    MOGA_Matsuoka_Run(whichCPG,whichGA_Case,trainingDataFile,[],i)
end

% clear all; 
% % close all; 
% clc;

%%


% % % % % general Matsuoka CPG case study:
whichCPG = '2_level_CPG3'; % CPG case study name
trainingDataFile = 'MatsRandomRes_6N_general_TrainingSet.mat'; % training data file name, for rea's reseacrch, hybrid GA/NN
generate_GenomeFile('2_level_CPG3'); % this functions holds definitions for the genetic sequence

%  GA only (With Feedback)
whichGA_Case = '_GA_only_Feedback';
%%
%  Begin Simulations:
for i=1:10  % running every method several times for stastistical analysis
    MOGA_Matsuoka_Run(whichCPG,whichGA_Case,trainingDataFile,[],i)
end