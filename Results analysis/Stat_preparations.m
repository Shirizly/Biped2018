addpath(genpath('results 0921'));
clear all
%%

len = [1:10];
mats = cell(1,length(len));
for i = 1:length(len)
%     mats{i} = ['6N_tagaLike_2Ank_torques_symm_feedback_eq_09_21_' num2str(len(i)) '.mat'];
    mats{i} = ['ConSpitz_eq_09_21_' num2str(len(i)) '.mat'];
%     mats{i} = ['2_level_CPG_eq_09_21_' num2str(len(i)) '.mat'];
end
if contains(mats{1},'tagaLike')
    casenum = 1;
else
    if contains(mats{1},'Spitz')
        casenum = 2;
    else
        if contains(mats{1},'2_level')
            casenum = 3;
        end
    end
end
    


stt = 2;
%%
for i = 1:numel(mats) %per Run
    clear GA
    load(mats{i});
    if GA.Progress<GA.Generations
        error(['bad file - ' num2str(i)]);
    end
    childstart = GA.Fittest(1)+GA.Fittest(2)+1;
    childend = GA.Population-GA.Fittest(3);
    childnum = childend-childstart+1;
    for gen = 1:GA.Generations % Per generation
        seqs = GA.Seqs(:,:,gen);
        Fits = GA.Fit(:,GA.FitIDs,gen);
        fronts = GA.get_fronts(GA,Fits);

        fitI = Fits(fronts{1},:);
        for f = 1:length(GA.FitIDs)
            Experts{i}(gen,f) = max(fitI(:,f));
            Averages{i}(gen,f) = mean(fitI(:,f));
        end
        %% Finding percentages of populations types not feasible
        fcount{i}(gen,1) = length(find(GA.Fit(:,1,gen)<1));
        fcountchild{i}(gen,1) = length(find(GA.Fit(childstart:childend,1,gen)<1));
        fcountrandom{i}(gen,1) = length(find(GA.Fit(childend+1:end,1,gen)<1));

        
        Afronts = GA.get_fronts(GA,Fits(:,2:3));
        AfitI = Fits(Afronts{1},2:3);
        scale = [1,12];
        
        Areas{i}(gen,1) = AreaUnderFront(AfitI,scale);
        
        %% low ankle torque:
        [TopIDs,~] = GA.GetTopPop(GA.Fittest(1));
        switch casenum
            case 1
                ind = find(seqs(TopIDs,4,end)./seqs(TopIDs,5,end)<0.05);
                lowAnkPer(i,gen) = length(ind)/GA.Fittest(1);
            case 2
                ind = find(abs(seqs(TopIDs,2,end)./max(seqs(TopIDs,[5,8],end),[],2))<0.05);
                lowAnkPer(i,gen) = length(ind)/GA.Fittest(1);
            case 3
                ind = find(abs(seqs(TopIDs,5,end)./seqs(TopIDs,8,end))<0.05);
                lowAnkPer(i,gen) = length(ind)/GA.Fittest(1);
        end
    end
end
%%
% figure
% plot(1:gen,cell2mat(Volumes))

%%
a_file = ['Stat_09_21_2' num2str(casenum)];
switch casenum
    case 1
        CPGName = 'TagaLike';
        CPGType = '6N_tagaLike_2Ank_torques_symm_feedback_eq'; 
    case 2
        CPGName = 'Spitz';
        CPGType = 'ConSpitz_eq';
    case 3  
        CPGName = 'Hybrid';
        CPGType = '2_level_CPG_eq';
end

a_full = [a_file '.mat'];
if exist(a_full,'file')==2 
    error('File already exists!')
else
    save(a_file,'Experts','Averages','CPGType','CPGName','fcount','fcountchild','fcountrandom','Areas','lowAnkPer');
    disp('file saved!');
end

%%
% col = [(0.1:0.2:1).',zeros(5,1),(1:-0.2:0.1).'];
% col2 = [(1:-0.2:0.1).',(0.1:0.2:1).',zeros(5,1)];
% col = [col;col2]
% col = [zeros(10,2),ones(10,1)];
col = [zeros(10,3)];
% col = [ones(10,1),zeros(10,2)];
figure
hold on
for i = 1:10
plot(1:20,lowAnkPer(i,:),'color',col(i,:),'lineWidth',1)
end
title(['Percentage of top-pop. with negligible ankle torque, CPG - '  CPGName]);
xlabel('Generation');
ylabel('Percentage');
axis([0 20 0 1])
grid on
hold off
