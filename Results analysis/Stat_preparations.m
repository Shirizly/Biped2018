addpath(genpath('results 0816'));
clear all
%%

len = [1:10];
mats = cell(1,length(len));
for i = 1:length(len)
    mats{i} = ['6N_tagaLike_2Ank_torques_symm_feedback_eq_08_14_' num2str(len(i)) '.mat'];
end



stt = 2;
%%
for i = 1:numel(mats) %per Run
    clear GA
    load(mats{i});
    childstart = GA.Fittest(1)+GA.Fittest(2)+1;
    childend = GA.Population-GA.Fittest(3);
    childnum = childend-childstart+1;
    for gen = 1:GA.Generations % Per generation
      
        Fits = GA.Fit(:,GA.FitIDs,gen);
        fronts = GA.get_fronts(GA,Fits);

        fitI = Fits(fronts{1},:);
        for f = 1:length(GA.FitIDs)
            Experts{i}(gen,f) = max(fitI(:,f));
            Averages{i}(gen,f) = mean(fitI(:,f));
        end
        
        fcount{i}(gen,1) = length(find(GA.Fit(:,1,gen)<1));
        fcountchild{i}(gen,1) = length(find(GA.Fit(childstart:childend,1,gen)<1));
        fcountrandom{i}(gen,1) = length(find(GA.Fit(childend+1:end,1,gen)<1));
        
        Volumes{i}(gen,1) = VolUnderFront(fitI);
    end
end
%%
% figure
% plot(1:gen,cell2mat(Volumes))

%%
a_file = 'Stat_08_16_02';
CPGName = 'TagaLike - Eq';
CPGType = '6N_tagaLike_2Ank_torques_symm_feedback_eq';
a_full = [a_file '.mat'];
if exist(a_full,'file')==2 
    error('File already exists!')
else
    save(a_file,'Experts','Averages','CPGType','CPGName','Volumes','fcount','fcountchild','fcountrandom');
    disp('file saved!');
end