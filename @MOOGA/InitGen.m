function [ GA ] = InitGen( GA )
%INIT Initializes the optimization variables
%   If an input MOOGA isn't provided, a random first generation is created

% Initialize sequences and fitness
GA.Seqs = zeros(GA.Population, GA.Gen.Length, GA.Generations);
GA.Parents = zeros(GA.Population, 2, GA.Generations); % 2 parents
GA.Fit = zeros(GA.Population, max(cell2mat(GA.FitFcn(:,1)')), GA.Generations);
GA.MLseqRunTime = NaN(GA.Population, 1, GA.Generations);
GA.simRunTime = NaN(GA.Population, 1, GA.Generations);
GA.sim_endCond = zeros(GA.Population, 1, GA.Generations);
GA.Tend_ratio = zeros(GA.Population, 1, GA.Generations);



GA.totGenTime = zeros(1,GA.Generations);
% TODO: do for 'MLseqRunTime' and 'simRunTime' what is done for 'Fit' and
%       'Seqs' (' if exist('In','var') == 1 ')

if exist(GA.FileIn,'file') == 2
    % Load input file
    In = load(GA.FileIn);
else
    if exist(GA.FileOut,'file') == 2
        % Load input file
        In = load(GA.FileOut);
    end
end

% if exist('In','var') == 1
%     % Are the same fitness functions used?
%     InNames = cell(In.GA.NFit,1);
%     ThisFitID = [];
%     InFitID = [];
%     for i = 1:In.GA.NFit
%         InNames{i} = MOOGA.GetFitFcnName(In.GA.FitFcn{i,2});
%     end
%     for i = 1:GA.NFit
%         ThisName = MOOGA.GetFitFcnName(GA.FitFcn{i,2});
%         ID = find(strcmp(ThisName,InNames),1,'first');
%         if ~isempty(ID)
%             ThisFitID = [ThisFitID,GA.FitFcn{i,1}]; %#ok<AGROW>
%             InFitID = [InFitID,In.GA.FitFcn{ID,1}]; %#ok<AGROW>
%         end
%         
%         if length(ThisFitID) ~= length(InFitID)
%             error('Discrepancy in FitFcn indexes')
%         end
%     end
% 
%     if In.GA.Population == GA.Population
%         if GA.ReDo
%             % Copy the last generation's seq. into new GA
%             % ids = [1,3:19];
%             % GA.Seqs(:,:,1) = In.GA.Seqs(:,ids,In.GA.Progress);
%             GA.Seqs(:,:,1) = In.GA.Seqs(:,:,In.GA.Progress);
%             GA.Fit(:,ThisFitID,1) = ...
%                 In.GA.Fit(:,InFitID,In.GA.Progress);
%             GA.Parents(:,:,1) = In.GA.Parents(:,:,In.GA.Progress);
%             GA.MLseqRunTime(:,1,1) = ...
%                 In.GA.MLseqRunTime(:,1,In.GA.Progress);
%             GA.simRunTime(:,1,1) = ...
%                 In.GA.simRunTime(:,1,In.GA.Progress);
%             GA.sim_endCond(:,1,1) = ...
%                 In.GA.sim_endCond(:,1,In.GA.Progress);
%             GA.Tend_ratio(:,1,1) = ...
%                 In.GA.Tend_ratio(:,1,In.GA.Progress);
%             GA.totGenTime(1,1) = ...
%                 In.GA.totGenTime(1,In.GA.Progress);
%                 
%         else
%             % Copy all the progress into new GA
%             GA.Seqs(:,:,1:In.GA.Progress) = ...
%                 In.GA.Seqs(:,:,1:In.GA.Progress);
%             GA.Fit(:,ThisFitID,1:In.GA.Progress) = ...
%                 In.GA.Fit(:,InFitID,1:In.GA.Progress);
%             GA.Parents(:,:,1) = ...
%                 In.GA.Parents(:,:,1:In.GA.Progress);
%             GA.MLseqRunTime(:,1,1) = ...
%                 In.GA.MLseqRunTime(:,1,1:In.GA.Progress);
%             GA.simRunTime(:,1,1) = ...
%                 In.GA.simRunTime(:,1,1:In.GA.Progress);
%             GA.sim_endCond(:,1,1) = ...
%                 In.GA.sim_endCond(:,1,1:In.GA.Progress);
%             GA.Tend_ratio(:,1,1) = ...
%                 In.GA.Tend_ratio(:,1,1:In.GA.Progress);
%             GA.totGenTime(1,1) = ...
%                 In.GA.totGenTime(1,1:In.GA.Progress);
%             
%             GA.Progress = In.GA.Progress-1;
%         end
%             
%     else
%         if In.GA.Population > GA.Population
%             warning('this code is not changed yet to do use bigger Pop!');
%             error('change it!');
%             
%             % Select the best genomes from the last generation
%             In.GA.Quant = GA.Quant;
%             TopIDs = In.GA.GetTopPop(GA.Population); % fitness = genes
% 
%             % Transfer top IDs to new population
%             GA.Seqs(:,:,1) = In.GA.Seqs(TopIDs,:,In.GA.Progress);  
%             if ~GA.ReDo 
%                 GA.Fit(:,ThisFitID,1) = ...
%                     In.GA.Fit(TopIDs,InFitID,In.GA.Progress);
%             end
%         else
%             % Copy the last generation's seq. into new GA
%             GA.Seqs(1:In.GA.Population,:,1) = ...
%                 In.GA.Seqs(:,:,In.GA.Progress);
%             if ~GA.ReDo
%                 GA.Fit(1:In.GA.Population,ThisFitID,1) = ...
%                     In.GA.Fit(:,InFitID,In.GA.Progress);
%             end
%             % Generate new random sequences
%             GA.Seqs(In.GA.Population+1:end,:,1) = ...
%                 GA.Gen.RandSeq(GA.Population-In.GA.Population);
%         end
%     end
% else


    % Generate new random sequences
    GA.Seqs(:,:,1) = GA.Gen.RandSeq(GA.Population);

    
    % for debugging purposes, loading the generation with working controllers
%         for i=1:GA.Population 
%           GA.Seqs(i,:,1) =  [ 1.3033  -12.3137   -0.9833    6.5217   -1.1679    0.5656    0.5375    0.4754    0.8111    0.1213    0.0083    0.3011    0.1157];
% 
% %         GA.Seqs(i,:,1) = [  0.1532   16.7269    0.5897    1.6174   -3.7263   10.4547   12.8639   -9.6789   18.1890    0.2962    0.2167    0.0791    0.0169];
%         end

% GA.Seqs(1,:,1) =[0.0638    2.0569   19.2291    1.4723    1.2078    0.0288    0.6410   -2.6102   -0.6777   18.9337];

end
                       


