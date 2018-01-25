
clc; close all; clear all;

load('VGAM_6N_TagaLike_symm_11_30_17_29_same_tonicInputs_20Gen_500Genes__NNclassi_and_rescale.mat')

inputNames = {'tau','b','w1','w2','w3','w4'};

% graphType = 'grey scale abs';
graphType = 'red and blue';
%% Show weights from input to the 1st layer:
NN_weights_matrix_plot( GA.NN_classi,{'IW',1},inputNames,graphType);

%% Show weights from layer j to layer i:
fromLayer = 1;
ToLayer = 2;
NN_weights_matrix_plot( GA.NN_classi,{'LW',ToLayer,fromLayer},inputNames,graphType);

%% also Show bias of layer i:
NN_weights_matrix_plot( GA.NN_classi,{'IW',1,'b'},inputNames,graphType);

%%
clc; close all;
NN_draw_graph( GA.NN_classi,inputNames);
