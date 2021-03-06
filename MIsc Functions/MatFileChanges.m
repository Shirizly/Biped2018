mats = {'Stat_04_25_1.mat',...
    'Stat_04_25_2.mat',...
    'Stat_04_25_3.mat'};
CPGns = {'2-level-v1',...
    '2-level-v2',...
    '2-level-v3'};
CPGs = {'2_level_CPG',...
    '2_level_CPG2',...
    '2_level_CPG3'};
for i = 1:numel(mats)
    clear Experts Averages
    load(mats{i})
    CPGType = CPGs{i};
    CPGName = CPGns{i};
    save(mats{i},'Experts','Averages','CPGType','CPGName');
end