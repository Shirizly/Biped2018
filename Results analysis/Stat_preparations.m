addpath(genpath('results 0921'),genpath('results 1115'));
clear all
%%

len = [1:10];
mats = cell(1,length(len));
for i = 1:length(len)
    mats{i} = ['6N_tagaLike_2Ank_torques_symm_feedback_eq_09_21_' num2str(len(i)) '.mat'];
%     mats{i} = ['ConSpitz_eq_11_15_' num2str(len(i)) '.mat'];
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
    mutestart = GA.Fittest(1)+1;
    mutenum = GA.Fittest(2)
    childstart = GA.Fittest(1)+GA.Fittest(2)+1;
    childend = GA.Population-GA.Fittest(3);
    childnum = childend-childstart+1
    randnum = GA.Fittest(3)
    for gen = 1:GA.Generations % Per generation
        seqs = GA.Seqs(:,:,gen);
        Fits = GA.Fit(:,GA.FitIDs,gen);
        fronts = GA.get_fronts(GA,Fits);

        fitI = Fits(fronts{1},:);
        for f = 1:length(GA.FitIDs)
            Experts{i}(gen,f) = max(fitI(:,f));
            Averages{i}(gen,f) = mean(fitI(:,f));
        end
        %% Finding percentages of populations types feasible
        fcount{i}(gen,1) = length(find(GA.Fit(:,5,gen)>0))/GA.Population;
        fcountmute{i}(gen,1) = length(find(GA.Fit(mutestart:childstart-1,5,gen)>0))/mutenum;
        fcountchild{i}(gen,1) = length(find(GA.Fit(childstart:childend,5,gen)>0))/childnum;
        fcountrandom{i}(gen,1) = length(find(GA.Fit(childend+1:end,5,gen)>0))/randnum;

        %%
        Afronts = GA.get_fronts(GA,Fits(:,2:3));
        AfitI = Fits(Afronts{1},2:3);
        scale = [1,12];
        
        Area1{i}(gen,1) = AreaUnderFront(AfitI,scale);
        %%
        Afronts = GA.get_fronts(GA,Fits(:,[1,3]));
        AfitI = Fits(Afronts{1},[1,3]);
        scale = [1,12];
        Area2{i}(gen,1) = AreaUnderFront(AfitI,scale);
        
        %%
        Afronts = GA.get_fronts(GA,Fits(:,1:2));
        AfitI = Fits(Afronts{1},1:2);
        scale = [1,12];
        Area3{i}(gen,1) = AreaUnderFront(AfitI,scale);
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

Areas = {Area1,Area2,Area3};
%%
% figure
% plot(1:gen,cell2mat(Volumes))

%%
a_file = ['Stat_09_21_5' num2str(casenum)];
switch casenum
    case 1
        CPGName = '1-level CPG';
        CPGType = '6N_tagaLike_2Ank_torques_symm_feedback_eq'; 
    case 2
        CPGName = '2-level CPG';
        CPGType = 'ConSpitz_eq';
    case 3  
        CPGName = 'Hybrid';
        CPGType = '2_level_CPG_eq';
end

a_full = [a_file '.mat'];
if exist(a_full,'file')==2 
    error('File already exists!')
else
    save(a_file,'Experts','Averages','CPGType','CPGName','fcount','fcountmute','fcountchild','fcountrandom','Areas','lowAnkPer');
    disp('file saved!');
end

%%
% col = [(0.1:0.2:1).',zeros(5,1),(1:-0.2:0.1).'];
% col2 = [(1:-0.2:0.1).',(0.1:0.2:1).',zeros(5,1)];
% col = [col;col2]
% col = [zeros(10,2),ones(10,1)];
% col = [zeros(10,3)];
% % col = [ones(10,1),zeros(10,2)];
% figure
% hold on
% for i = 1:10
% plot(1:20,lowAnkPer(i,:),'color',col(i,:),'lineWidth',1)
% end
% title(['Percentage of top-pop. with negligible ankle torque, CPG - '  CPGName]);
% xlabel('Generation');
% ylabel('Percentage');
% axis([0 20 0 1])
% grid on
% hold off
