function [  ] =MOGA_NN_crossvalidation(whichCPG,whichNN,trainingDataFile)
% MOGA_MATSUOKA_RUN - run MOGA optimization
% INPUTS:
% *) 'whichCPG' - which CPG to run?
%                 #) '2N_symm' - 2-nueron symmetric matsuoka CPG
%                 #) '2N_general' - 2-nueron general matsuoka CPG
%                 #) '6N_TagaLike_general' - 6-neurons general coupling weights between the hip and the ankles
%                 #) '6N_TagaLike_symm' - 6-neurons symmetric coupling weights between the hip and the ankles
%                 #) '6N_TagaLike_with_hip_FB' 
% *) 'whichNN' - 'regression' or 'classification'
% *) 'trainingDataFile' - name of the NN training data file
% *) 'prevMOGAfileIn' - name of previuos MOGA run for continuing optimization
 
% get MOGA output file name:
FileName_prefix = 'NN_crossValidation_';
FileName_start = [whichCPG,'_'];
FileName_date = datestr(now,'mm_dd_hh_MM');
FileName_extra1 = '_10fold_';
FileOut_name = [FileName_prefix,FileName_start,FileName_date,...
    FileName_extra1,'.mat'];

% Load Genome File:
load('MatsuokaGenome.mat','Keys','Range','N',...
    'nAnkle1','nAnkle2','nHip','maxAnkle', 'maxHip',...
    'MutDelta0','MutDelta1');

Gen = Genome(Keys, Range);

MML = MatsuokaML();
% MML.perLim = [0.68 0.78];
% MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.perLim = [1.3 1.5];
MML.perLimOut = MML.perLim + [-0.08 0.08]; % Desired period range
MML.nNeurons = 2*N;
MML.normParams = [];

% % % define the CPG weights names:
switch 'whichCPG'
    case '2N_symm'
        W_name = '2N_symm'; %correct this
        neuronNumInCPG = '2N_CPG';
    case '2N_general'
        W_name = '2N_general'; %correct this
        neuronNumInCPG = '2N_CPG';
    case '6N_TagaLike_general'
        W_name = '6neuron_taga_like';
        neuronNumInCPG = '6N_CPG';
    case '6N_TagaLike_symm'
        W_name = '6neuron_taga_like_symm';
        neuronNumInCPG = '6N_CPG';
    case '6N_general'
        W_name = 'weights';
        neuronNumInCPG = '6N_CPG';
    otherwise
        error('unknown CPG type...');
end


maxN = 250000;
inFilenames = {trainingDataFile};
switch whichNN % preparing the inputs and outputs names
    case 'regression'
        MML.sample_genes = {'\tau_r','beta',W_name}; 
        MML.target_genes = {'period'};
    case 'classification'
        MML.sample_genes = {'\tau_r','beta',W_name}; 
        MML.target_genes = {'n_osc and osc classes'}; % not in use
end

% % run k-fold crossvalidation:
Kfold_num = 10;
NNtypes = {[5],[10],[15],[20],[25],[30],[35],[40],[45],[50]};
NN_perf_store_training = zeros(Kfold_num,numel(NNtypes));
NN_perf_store_validation = zeros(Kfold_num,numel(NNtypes));
NN_perf_store_testing = zeros(Kfold_num,numel(NNtypes));

for i=1:numel(NNtypes)
    for k=1:Kfold_num
        architecture = NNtypes{1,k};
        switch whichNN % prepare and train the neural network
            case 'regression'
                [samples, targets, ~] = MML.prepare_reg_NNData(neuronNumInCPG,inFilenames, maxN);

                [~, tr] = ...
                        MML.train_reg_NN(samples, targets, architecture, NNSamples);

            case 'classification'
                [samples, targets] = ...
                    MML.prepare_classi_NNData(neuronNumInCPG,inFilenames, maxN);

                [~, tr] = MML.train_classi_NN(samples, targets, architecture);
        end

        % % save the training performance:
        NN_perf_store_training(k,i) = tr.perf;
        NN_perf_store_validation(k,i) = tr.perf_val;
        NN_perf_store_testing(k,i) = tr.perf_test;
    end
end

save(FileOut_name,'Kfold_num','NNtypes','NN_perf_store_training',...
    'NN_perf_store_validation','NN_perf_store_testing');

% TODO: plotting functions

end