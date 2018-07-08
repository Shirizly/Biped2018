clear all



% mats = {'AS_2_level_CPG_06_06_1.mat';
%     'AS_2_level_CPG_06_06_2.mat';
%     'AS_2_level_CPG_06_06_3.mat';
%     'AS_2_level_CPG_06_06_4.mat';
%     'AS_2_level_CPG_06_06_5.mat';
%     'AS_2_level_CPG_06_06_6.mat';
%     'AS_2_level_CPG_06_06_7.mat';
%     'AS_2_level_CPG_06_06_8.mat';
%     'AS_2_level_CPG_06_06_9.mat';
%     'AS_2_level_CPG_06_06_10.mat'};

% mats = {'AS_2_level_CPG2_06_04_1.mat';
%     'AS_2_level_CPG2_06_04_2.mat';
%     'AS_2_level_CPG2_06_04_3.mat';
%     'AS_2_level_CPG2_06_05_4.mat';
%     'AS_2_level_CPG2_06_05_5.mat';
%     'AS_2_level_CPG2_06_05_6.mat';
%     'AS_2_level_CPG2_06_05_7.mat';
%     'AS_2_level_CPG2_06_05_8.mat';
%     'AS_2_level_CPG2_06_05_9.mat';
%     'AS_2_level_CPG2_06_05_10.mat'};

% mats = {'AS_2_level_CPG3_06_05_1.mat';
%     'AS_2_level_CPG3_06_05_2.mat';
%     'AS_2_level_CPG3_06_05_3.mat';
%     'AS_2_level_CPG3_06_05_4.mat';
%     'AS_2_level_CPG3_06_05_5.mat';
%     'AS_2_level_CPG3_06_05_6.mat';
%     'AS_2_level_CPG3_06_06_7.mat';
%     'AS_2_level_CPG3_06_06_8.mat';
%     'AS_2_level_CPG3_06_06_9.mat';
%     'AS_2_level_CPG3_06_06_10.mat'};

% mats = {'AS_6N_tagaLike_2Ank_torques_symm_feedback_06_05_1.mat',...
%     'AS_6N_tagaLike_2Ank_torques_symm_feedback_06_05_2.mat',...
%     'AS_6N_tagaLike_2Ank_torques_symm_feedback_06_05_3.mat',...
%     'AS_6N_tagaLike_2Ank_torques_symm_feedback_06_05_4.mat',...
%     'AS_6N_tagaLike_2Ank_torques_symm_feedback_06_06_5.mat',...
%     'AS_6N_tagaLike_2Ank_torques_symm_feedback_06_06_6.mat',...
%     'AS_6N_tagaLike_2Ank_torques_symm_feedback_06_06_7.mat',...
%     'AS_6N_tagaLike_2Ank_torques_symm_feedback_06_09_8.mat',...
%     'AS_6N_tagaLike_2Ank_torques_symm_feedback_06_10_9.mat',...
%     'AS_6N_tagaLike_2Ank_torques_symm_feedback_06_10_10.mat'};

% mats = {'AS_ConSpitz_06_15_1.mat',...
%     'AS_ConSpitz_06_15_2.mat',...
%     'AS_ConSpitz_06_15_3.mat',...
%     'AS_ConSpitz_06_15_4.mat',...
%     'AS_ConSpitz_06_15_5.mat',...
%     'AS_ConSpitz_06_15_6.mat',...
%     'AS_ConSpitz_06_15_7.mat',...
%     'AS_ConSpitz_06_15_8.mat',...
%     'AS_ConSpitz_06_15_9.mat',...
%     'AS_ConSpitz_06_15_10.mat'};

% mats = {'AS_ConSpitz2_06_15_1.mat',...
%     'AS_ConSpitz2_06_15_2.mat',...
%     'AS_ConSpitz2_06_15_3.mat',...
%     'AS_ConSpitz2_06_15_4.mat',...
%     'AS_ConSpitz2_06_15_5.mat',...
%     'AS_ConSpitz2_06_15_6.mat',...
%     'AS_ConSpitz2_06_15_7.mat',...
%     'AS_ConSpitz2_06_15_8.mat',...
%     'AS_ConSpitz2_06_15_9.mat',...
%     'AS_ConSpitz2_06_15_10.mat'};

mats = {'AS_ConSpitz3_06_14_1.mat',...
    'AS_ConSpitz3_06_14_2.mat',...
    'AS_ConSpitz3_06_14_3.mat',...
    'AS_ConSpitz3_06_14_4.mat',...
    'AS_ConSpitz3_06_15_5.mat',...
    'AS_ConSpitz3_06_15_6.mat',...
    'AS_ConSpitz3_06_15_7.mat',...
    'AS_ConSpitz3_06_15_8.mat',...
    'AS_ConSpitz3_06_15_9.mat',...
    'AS_ConSpitz3_06_15_10.mat'};
%%
for i = 1:numel(mats) %per Run
    clear GA
    load(mats{i});
    for gen = 1:GA.Generations % Per generation
      
        Fits = GA.Fit(:,GA.FitIDs,gen);
        fronts = GA.get_fronts(GA,Fits);

        fitI = Fits(fronts{1},:);
        for f = 1:length(GA.FitIDs)
            Experts{i}(gen,f) = max(fitI(:,f));
            Averages{i}(gen,f) = mean(fitI(:,f));
        end
    end
end


%%
a_file = 'Stat_06_13_7';
CPGName = 'Spitz - FAM2';
CPGType = 'ConSpitz3';
a_full = [a_file '.mat'];
if exist(a_full,'file')==2 
    error('File already exists!')
else
    save(a_file,'Experts','Averages','CPGType','CPGName');
    disp('file saved!');
end