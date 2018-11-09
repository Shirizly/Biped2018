function CollectTop(ColName,GANames,CritType,varargin)
% receives as input:
% collection file name
% cell array of names of GA results files
% Criterion type: 
% 1 - Experts
% 2 - 1st front
% 3 - Top pop
% 4 - JOAT (decide how!!)

if nargin>3
    Gens = varargin{1};
else
    Gens = 

if exist([ColName '.mat'],'file') == 2
    load([ColName '.mat']);
end

for i = 1:size(GANames,1)
    load([GANames{i} '.mat'])
    for 

end
