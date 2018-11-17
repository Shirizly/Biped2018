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
% addpath(genpath('results 0902'))


%%


flag = 1;

mats = {'Stat_09_21_21.mat',...
    'Stat_09_21_22.mat',...
    'Stat_09_21_23.mat'};

mats = {'Stat_09_21_21.mat',...
    'Stat_09_21_22.mat'};
n = max(size(mats,2),size(mats,1));
col = {'r','b','k','g','c'};
colind = [1,0,0;0,0,1;0,0,0];
alpgen = 0.05;

%%
alp = alpgen;
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


for j = size(ExpertVel{i},2):-1:1
    maxval = max(expmean(:,j));
    maxind = find(expmean(:,j)==maxval);
    range = [1:maxind-1 maxind+1:numel(mats)];
    secbestval = max(expmean(range,j));
    secbestind = find(expmean(:,j)==secbestval);
    [p(j),h(j)] = ranksum(ExpertVel{maxind}(:,j),ExpertVel{secbestind}(:,j),'tail','right','alpha',alp);
    colv = [1,1,1];
    if h(j) == 1
        colv = colind(maxind,:);
    end
    if colv == [1,1,1]
        pl(end+1) = plot(j,-0,'o','MarkerSize',10,'MarkerEdgeColor','k');
    else
        pl(end+1) = plot(j,-0,'o','MarkerFaceColor',colv,'MarkerSize',10,'MarkerEdgeColor','k');
    end
end



leg{n+1}= ['RankSumTest - ' num2str((1-alp)*100) '%'];


legend(leg,'location','best','fontSize',12,'FontWeight','bold')
xlabel('Generation','FontName','Garamond','FontWeight','bold','fontSize',18)
ylabel('Velocity FF','FontName','Garamond','FontWeight','bold','fontSize',18)
title('Velocity Experts GA progression','FontName','Garamond','FontWeight','bold','fontSize',22)
grid on


hold off




%%
alp = alpgen;
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
    expmean(i,:) = y;
%     bootfun = @(x)(median(x));

    erbr = bootci(1000,{bootfun,ExpertSto{i}},'alpha',alp);
    pl(i) = errorbar(x,y,y-erbr(1,:),erbr(2,:)-y,'Color',col{i},'LineWidth',2);
    leg{i} = CPGName;
    
end

for j = size(ExpertSto{1},2):-1:1
    maxval = max(expmean(:,j));
    maxind = find(expmean(:,j)==maxval);
    range = [1:maxind-1 maxind+1:numel(mats)];
    secbestval = max(expmean(range,j));
    secbestind = find(expmean(:,j)==secbestval);
    [p(j),h(j)] = ranksum(ExpertSto{maxind}(:,j),ExpertSto{secbestind}(:,j),'tail','right','alpha',alp);
    colv = [1,1,1];
    if h(j) == 1
        colv = colind(maxind,:);
    end
    if colv == [1,1,1]
        pl(end+1) = plot(j,-0,'o','MarkerSize',10,'MarkerEdgeColor','k');
    else
        pl(end+1) = plot(j,-0,'o','MarkerFaceColor',colv,'MarkerSize',10,'MarkerEdgeColor','k');
    end

end





leg{n+1}= ['SignRankTest - ' num2str((1-alp)*100) '%'];
axis([0 20 0 1])
legend(leg,'location','best','fontSize',12,'FontWeight','bold')
xlabel('Generation','FontName','Garamond','FontWeight','bold','fontSize',18)
ylabel('Stochastic FF','FontName','Garamond','FontWeight','bold','fontSize',18)
title('Rough Terrain Experts GA progression','FontName','Garamond','FontWeight','bold','fontSize',22)
grid on


hold off


%%
alp = alpgen;
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
    expmean(i,:) = y;
%     bootfun = @(x)(median(x));

    erbr = bootci(1000,{bootfun,ExpertUDH{i}},'alpha',0.05);
    pl(i) = errorbar(x,y,y-erbr(1,:),erbr(2,:)-y,'Color',col{i},'LineWidth',2);
    leg{i} = CPGName;
    
end

for j = 1:size(ExpertUDH{i},2)
    maxval = max(expmean(:,j));
    maxind = find(expmean(:,j)==maxval);
    range = [1:maxind-1 maxind+1:numel(mats)];
    secbestval = max(expmean(range,j));
    secbestind = find(expmean(:,j)==secbestval);
    [p(j),h(j)] = ranksum(ExpertUDH{maxind}(:,j),ExpertUDH{secbestind}(:,j),'tail','right','alpha',0.01);
    colv = [1,1,1];
    if h(j) == 1
        colv = colind(maxind,:);
    end
    if colv == [1,1,1]
        pl(end+1) = plot(j,-0,'o','MarkerSize',10,'MarkerEdgeColor','k');
    else
        pl(end+1) = plot(j,-0,'o','MarkerFaceColor',colv,'MarkerSize',10,'MarkerEdgeColor','k');
    end
end

leg{n+1}= ['RankSumTest - ' num2str((1-alp)*100) '%'];


legend(leg,'location','best','fontSize',12,'FontWeight','bold')
xlabel('Generation','FontName','Garamond','FontWeight','bold','fontSize',18)
ylabel('Slope Robustness FF','FontName','Garamond','FontWeight','bold','fontSize',18)
title('Slope Robustness Experts GA progression','FontName','Garamond','FontWeight','bold','fontSize',22)

grid on
hold off

%%
alp = alpgen;
leg = [];
figure
hold on
for i = 1:numel(mats)
    load(mats{i})
    for j=1:numel(Areas)
        Ar{i}(j,:) = Areas{j}(:).';

    end
    x = 1:size(Ar{i},2);
    bootfun = @(x)(mean(x));
    y = median(Ar{i});
    Armean(i,:) = y;
%     bootfun = @(x)(median(x));

    erbr = bootci(1000,{bootfun,Ar{i}},'alpha',0.05);
    pl(i) = errorbar(x,y,y-erbr(1,:),erbr(2,:)-y,'Color',col{i},'LineWidth',2);
    leg{i} = CPGName;
    
end

for j = 1:size(Ar{i},2)
    if all(Armean(:,j)==0)
        pl(end+1) = plot(j,-0,'o','MarkerSize',10,'MarkerEdgeColor','k');
    else
        maxval = max(Armean(:,j));
        maxind = find(Armean(:,j)==maxval);
        range = [1:maxind-1 maxind+1:numel(mats)];
        secbestval = max(Armean(range,j));
        secbestind = find(Armean(:,j)==secbestval);
        [p(j),h(j)] = signrank(Ar{maxind}(:,j),Ar{secbestind}(:,j),'tail','right','alpha',0.05);
        colv = [1,1,1];
        if h(j) == 1
            colv = colind(maxind,:);
        end
        if colv == [1,1,1]
            pl(end+1) = plot(j,-0,'o','MarkerSize',10,'MarkerEdgeColor','k');
        else
            pl(end+1) = plot(j,-0,'o','MarkerFaceColor',colv,'MarkerSize',10,'MarkerEdgeColor','k');
        end
    end
end

leg{n+1}= ['SignRankTest - ' num2str((1-alp)*100) '%'];


legend(leg,'location','best','fontSize',12,'FontWeight','bold')
xlabel('Generation','fontSize',18,'FontName','Garamond','FontWeight','bold')
ylabel('Area under 1st Front','FontName','Garamond','fontSize',18,'FontWeight','bold')
title('Area under 1st Front GA progression','FontName','Garamond','fontSize',22,'FontWeight','bold')

grid on
hold off