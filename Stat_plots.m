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



%% preparing the data into vectors for  plots
% for i=1:numel(Experts)
% ExpertVel(i,:) = Experts{i}(:,1).';
% ExpertEng(i,:) = Experts{i}(:,2).';
% end

%%


flag = 1;
mats = {'Stat_06_13_4.mat',...
        'Stat_06_13_5.mat',...
    'Stat_06_13_1.mat'};

% mats = {'Stat_06_13_5.mat',...
%     'Stat_06_13_6.mat',...
%     'Stat_06_13_7.mat'};
n = size(mats,2);
col = {'b','r','k','g','c'};

%%
% hold on
% for i = 1:numel(mats)
%     load(mats{i})
%     for j=1:numel(Experts)
%         ExpertVel{i}(j,:) = Experts{j}(:,1).';
%         ExpertVel{i}(j,:) = max(ExpertVel{i}(j,:)-1,0)
% %         ExpertSto(i,:) = Experts{i}(:,2).';
%     end
%     x = 1:size(ExpertVel{i},2);
%     bootfun = @(x)(mean(x));
%     y = median(ExpertVel{i});
%     
% %     bootfun = @(x)(median(x));
% 
%     erbr = bootci(1000,{bootfun,ExpertVel{i}},'alpha',0.05);
%     pl(i) = errorbar(x,y,y-erbr(1,:),erbr(2,:)-y,'Color',col{i},'LineWidth',2);
%     leg{i} = CPGName;
%     
% end
% leg
% 
% for j = 1:size(ExpertVel{i},2)
% [p(j),h(j)] = ranksum(ExpertVel{5}(:,j),ExpertVel{1}(:,j));
% col1 = [1-h(j),0,h(j)];
% pl(end+1) = plot(j,0,'o','MarkerFaceColor',col1,'MarkerSize',10,'MarkerEdgeColor','k');
% end
% 
% % plot(1:size(ExpertVel{i},2),h,'*')
% leg{n+1}= 'RankSumTest-v5/v1';
% legend(leg)
% xlabel('Generation')
% ylabel('Velocity FF')
% title('Velocity Experts GA progression - over 5 runs')
% grid on
% 
% 
% hold off

%%
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

for j = size(ExpertSto{i},2):-1:1
[p(j),h(j)] = ranksum(ExpertSto{1}(:,j),ExpertSto{2}(:,j),'tail','left','method','exact');
col = [1-h(j),0,h(j)];
pl(end+1) = plot(j,0,'o','MarkerFaceColor',col,'MarkerSize',10,'MarkerEdgeColor','k');



[p1(j),h1(j)] = ranksum(ExpertSto{3}(:,j),ExpertSto{2}(:,j),'tail','left','method','exact');
col = [1-h1(j),0,h1(j)];
pl(end+1) = plot(j,0.1,'gs','MarkerFaceColor',col,'MarkerSize',10,'MarkerEdgeColor','k');
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
%%
% 
% for i = 1:numel(mats)
%     load(mats{i})
%     for j=1:numel(Experts)
%         ExpertSto{i}(j,:) = Experts{j}(:,2).';
% 
%     end
%     x = 1:size(ExpertSto{i},2);
%     bootfun = @(x)(mean(x));
%     y = median(ExpertSto{i});
%     
% %     bootfun = @(x)(median(x));
% 
%     erbr = bootci(1000,{bootfun,ExpertSto{i}},'alpha',0.05);
%     pl(i) = errorbar(x,y,y-erbr(1,:),erbr(2,:)-y,'Color',col{i},'LineWidth',2);
%     leg{i} = CPGName;
%     
% end
% 
% for j = 1:size(ExpertSto{i},2)
% [p(j),h(j)] = ranksum(ExpertSto{2}(:,j),ExpertSto{1}(:,j));
% col = [1-h(j),0,h(j)];
% pl(end+1) = plot(j,0.8,'o','MarkerFaceColor',col,'MarkerSize',10,'MarkerEdgeColor','k');
% end
% 
% leg{4}= 'RankSumTest-v2/v1';
% legend(leg)
% xlabel('Generation')
% ylabel('Energy efficiency FF')
% title('Energy efficiency Experts GA progression - over 5 runs')
% axis([0 20 0.8 1])
% grid on
% 
% 
% hold off