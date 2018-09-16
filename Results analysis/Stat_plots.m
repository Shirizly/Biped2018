%% simple plots of the data
% figure
% hold on
% for i=1:10
% plot((1:size(Experts{i},1)),Experts{i}(:,1))
% end
% hold off
% figure
% hold on
% for i=1:10
% plot((1:size(Experts{i},1)),Experts{i}(:,2))
% end
% hold off
%%
clc
clear



%%


flag = 1;
mats = {'Stat_08_16_01.mat',...
    'Stat_08_16_02.mat',...
    'Stat_08_16_03.mat'};

% mats = {'Stat_09_02_1.mat','Stat_08_16_1.mat'};
% mats = {'Stat_06_13_5.mat',...
%     'Stat_06_13_6.mat',...
%     'Stat_06_13_7.mat'};
n = size(mats,2);
col = {'b','r','k','g','c'};

%% Checking for zero values (GA didn't finish)
% load(mats{3})
% figure
% hold on
% for j=1:numel(Experts)
%     ExpertVel(j,:) = Experts{j}(:,1).';
% %     ExpertVel(j,:) = max(ExpertVel(j,:)-1,0)
%     
%     
% end
% for j = 1:size(ExpertVel,2)
%     plot(j*ones(size(ExpertVel(:,j))),ExpertVel(:,j))
% end
% hold off
%%
leg = [];
figure
hold on
for i = 1:numel(mats)
    load(mats{i})
    for j=1:numel(Experts)
        ExpertVel{i}(j,:) = Experts{j}(:,1).';
        ExpertVel{i}(j,:) = max(ExpertVel{i}(j,:)-1,0);

    end
    x = 1:size(ExpertVel{i},2);
    bootfun = @(x)(mean(x));
    y = mean(ExpertVel{i});
    expmean(i,:) = y;
%     bootfun = @(x)(median(x));

    erbr = bootci(1000,{bootfun,ExpertVel{i}},'alpha',0.05);
    pl(i) = errorbar(x,y,y-erbr(1,:),erbr(2,:)-y,'Color',col{i},'LineWidth',2);
    leg{i} = CPGName;
    
end
% leg

for j = size(ExpertVel{i},2):-1:1
[p(j),h(j)] = ranksum(ExpertVel{1}(:,j),ExpertVel{2}(:,j),'tail','right');
col1 = [1-h(j),0,h(j)];
pl(end+1) = plot(j,-0,'o','MarkerFaceColor',col1,'MarkerSize',10,'MarkerEdgeColor','k');

% [p1(j),h1(j)] = ranksum(ExpertVel{1}(:,j),ExpertVel{3}(:,j),'tail','right');
% col2 = [1-h1(j),0,h1(j)];
% pl(end+1) = plot(j,-0.05,'gs','MarkerFaceColor',col2,'MarkerSize',10,'MarkerEdgeColor','k');
end
% [pr(j),hr(j)] = signrank(expmean(3,:),expmean(2,:),'tail','right');
% colr = [1-hr(j),0,hr(j)];
% pl(end+1) = plot(0,0.5,'gs','MarkerFaceColor',colr,'MarkerSize',15,'MarkerEdgeColor','k');


% plot(1:size(ExpertVel{i},2),h,'*')
leg{n+1}= 'RankSumTest-S1/S2';
% leg{n+2}= 'SignRankTest-S1/S2';
legend(leg)
xlabel('Generation')
ylabel('Velocity FF')
title('Velocity Experts GA progression - over 10 runs')
grid on


hold off




%%
leg = [];
figure
hold on
for i = 1:numel(mats)
    load(mats{i})
    for j=1:numel(Experts)
        ExpertSto{i}(j,:) = Experts{j}(:,2).';

    end
    x = 1:size(ExpertSto{i},2);
    bootfun = @(x)(mean(x));
    y = median(ExpertSto{i});
    
%     bootfun = @(x)(median(x));

    erbr = bootci(1000,{bootfun,ExpertSto{i}},'alpha',0.05);
    pl(i) = errorbar(x,y,y-erbr(1,:),erbr(2,:)-y,'Color',col{i},'LineWidth',2);
    leg{i} = CPGName;
    
end

for j = size(ExpertSto{1},2):-1:1
    comp1 = ExpertSto{1}(:,j);
    comp2 = ExpertSto{2}(:,j);
    comp3 = ExpertSto{3}(:,j);
    
    erbr = bootci(1000,{bootfun,comp1},'alpha',0.05);
    remove1 = [find(comp1>erbr(2));0;find(comp1<erbr(1))].';
%     comp1(remove1) = [];
    erbr = bootci(1000,{bootfun,comp2},'alpha',0.05);
    remove2 = [find(comp2>erbr(2));find(comp2<erbr(1))].';
%     
% 
%     simel = [];
%     for c1 = 1:length(remove1)
%         findel = find(remove2 == remove1(c1));
%         if ~isempty(findel)
%             simel = [simel remove2(findel)];
%         end
%     end
%     if ~isempty(simel)
%         disp(simel)
%         comp1(simel) = [];
%         comp2(simel) = [];
%     end
    
    dev(j) = std(comp1);
    dev2(j) = std(comp2);
    var3 = std(comp3);
    mean1 = mean(comp1);
    mean2 = mean(comp2);
    mean3 = mean(comp3);
    
%     disp([dev(j);dev2(j);dev3(j)])
    
%     [p(j),h(j)] = ranksum(comp3,comp2,'tail','right','method','exact');
    [p1(j),h1(j)] = signrank(comp2,comp1,'tail','right','method','exact');
    colb = [1-h(j),0,h(j)];
% 
    pl(end+1) = plot(j,-0.1,'o','MarkerFaceColor',colb,'MarkerSize',10,'MarkerEdgeColor','k');
%     
%     pl(end+1) = plot(j,p(j),'o','MarkerFaceColor',colb,'MarkerSize',10,'MarkerEdgeColor','k');
% % %     
    [p1(j),h1(j)] = signrank(comp2,comp3,'tail','right','method','exact');
    colc = [1-h1(j),0,h1(j)];
    pl(end+1) = plot(j,-0.05,'gs','MarkerFaceColor',colc,'MarkerSize',10,'MarkerEdgeColor','k');
%     pl(end+1) = plot(j,p1(j),'gs','MarkerSize',10,'MarkerFaceColor',colc);

end


leg{n+1}= 'RankSumTest-S2/S1';
leg{n+2}= 'RankSumTest-S2/S3';
legend(leg)
xlabel('Generation')
ylabel('Stochastic FF')
title('Stochastic Experts GA progression - over 10 runs')
grid on


hold off


%%
leg = [];
figure
hold on
for i = 1:numel(mats)
    load(mats{i})
    for j=1:numel(Experts)
        ExpertUDH{i}(j,:) = Experts{j}(:,3).';

    end
    x = 1:size(ExpertUDH{i},2);
    bootfun = @(x)(mean(x));
    y = median(ExpertUDH{i});
    
%     bootfun = @(x)(median(x));

    erbr = bootci(1000,{bootfun,ExpertUDH{i}},'alpha',0.05);
    pl(i) = errorbar(x,y,y-erbr(1,:),erbr(2,:)-y,'Color',col{i},'LineWidth',2);
    leg{i} = CPGName;
    
end

for j = 1:size(ExpertUDH{i},2)
[p(j),h(j)] = ranksum(ExpertUDH{1}(:,j),ExpertUDH{2}(:,j),'tail','right');
% [p(j),h(j)] = ranksum(ExpertUDH{1}(:,j),ExpertUDH{2}(:,j),'tail','right');
col2 = [1-h(j),0,h(j)];
pl(end+1) = plot(j,p(j),'o','MarkerFaceColor',col2,'MarkerSize',10,'MarkerEdgeColor','k');
end

leg{n+1}= 'RankSumTest-S1/S2';


legend(leg)
xlabel('Generation')
ylabel('Slope Robustness FF')
title('Slope Robustness Experts GA progression - over 10 runs')

grid on
hold off

%%
leg = [];
figure
hold on
for i = 1:numel(mats)
    load(mats{i})
    for j=1:numel(Volumes)
        Vol{i}(j,:) = Volumes{j}(:).';

    end
    x = 1:size(Vol{i},2);
    bootfun = @(x)(mean(x));
    y = median(Vol{i});
    
%     bootfun = @(x)(median(x));

    erbr = bootci(1000,{bootfun,Vol{i}},'alpha',0.05);
    pl(i) = errorbar(x,y,y-erbr(1,:),erbr(2,:)-y,'Color',col{i},'LineWidth',2);
    leg{i} = CPGName;
    
end

for j = 1:size(Vol{i},2)
[p(j),h(j)] = signrank(Vol{2}(:,j),Vol{1}(:,j),'tail','right');
% [p(j),h(j)] = signrank(ExpertUDH{1}(:,j),ExpertUDH{2}(:,j),'tail','right');
col2 = [1-h(j),0,h(j)];
pl(end+1) = plot(j,-0.05,'o','MarkerFaceColor',col2,'MarkerSize',10,'MarkerEdgeColor','k');
end

leg{n+1}= 'RankSumTest-S2/S1';


legend(leg)
xlabel('Generation')
ylabel('Volume under 1st Front')
title('Volume under 1st Front GA progression - over 10 runs')

grid on
hold off