function [  ] =...
    MOGA_Matsuoka_Run(whichCPG,whichGA_Case,trainingDataFile,prevMOGAfileIn,startdate,runcount)
% MOGA_MATSUOKA_RUN - run MOGA optimization
% INPUTS:
% *) 'whichCPG' - which CPG to run?
%                 #) '2N_symm' - 2-nueron symmetric matsuoka CPG
%                 #) '2N_general' - 2-nueron general matsuoka CPG
%                 #) '6N_TagaLike_general' - 6-neurons general coupling weights between the hip and the ankles
%                 #) '6N_TagaLike_symm' - 6-neurons symmetric coupling weights between the hip and the ankles
%                 #) '6N_TagaLike_with_hip_FB' 
% *) 'whichGA_Case' - with or without 'NN', rescale ect.
%                       '_GA_only','_rescale_only','_NN_classi_only',
% *) 'trainingDataFile' - name of the NN training data file
% *) 'prevMOGAfileIn' - name of previuos MOGA run for continuing optimization
                
% Evulotionary parameters:
MoogaGen = 20;
MoogaPop = 500;
GA = MOOGA(MoogaGen,MoogaPop); % (generations,population size)
GA = GA.SetFittest(15,15,2); %(%keep, %mutate, %random), all the rest are created with crossing
% GA = GA.SetFittest(100,0,0);
GA.JOAT = 0; % jack of all trades, methods to prevent overspecializations for specific fitnesses (see jonathan's thesis)
GA.Quant = 0.7; % parameter for GA.JOAT = 2, defines percentiles for JOAT method 2

GA.FileIn = prevMOGAfileIn; % load previous mooga result file, can continue from stop point

% get MOGA output file name:
FileName_prefix = 'AS_';
FileName_start = [whichCPG,'_'];
FileName_date = [startdate '_' num2str(runcount)];
FileName_extra1 = [];
FileName_extra2 = [];
% FileName_extra2 = whichGA_Case;
GA.FileOut = [FileName_prefix,FileName_start,FileName_date,...
    FileName_extra1,FileName_extra2,'.mat'];

% Load Genome File:
load('MatsuokaGenome.mat','Keys','Range','N',...
    'nAnkle1','nAnkle2','nHip','maxAnkle', 'maxHip',...
    'MutDelta0','MutDelta1');



switch whichGA_Case
    case {'_GA_only','_GA_only_Feedback'}
        % Use NN?
        use_NN = 0;

   case '_rescale_only'
        use_NN = 0;
        GA.rescaleFcn = @rescaleFcn; 
        
    case {'_NN_classi_only','_NN_classi_only_Feedback'}'
        use_NN = 'NN_classi';
    
    case '_NNs_and_rescale'
        % TODO: incoporate both REG and CLASSI NNs in the code
        use_NN = 'NN_reg';
        GA.rescaleFcn = @rescaleFcn;
        
    case '_NNclassi_and_rescale'
        use_NN = 'NN_classi';
        GA.rescaleFcn = @rescaleFcn;
end
             
GA.Graphics = 0;
GA.ReDo = 1;

GA.Gen = Genome(Keys, Range);



%%
% function noGO_flag = genomeChevkFcn(X,T)
%     % If the ankle torque is a posivite constant (larger than zero and with
%     % no period) than don't use this genome.
%     [~, periods, signals, ~, ~] = MML.processResults(X, T);
%     
%     % check that the ankle torque is positive:
%     ind2analize = ceil(0.7*size(signals,2));
%     lastOfAnkleSignal = signals(1,ind2analize:end);
%     
%     %  dont take CPGs if the ankle torque is not oscillating and positive
%     if (mean(lastOfAnkleSignal) > 0.001) && isnan(periods(1,1))
%         noGO_flag = 1;
%     else % if the hip torque is constant than ignore it as well
%         if isnan(periods(2,1))
%             noGO_flag = 1;
%         else
%             noGO_flag = 0;
%         end
%     end


%% Set up the simulations
GA.Sim = Simulation();
GA.Sim.Graphics = GA.Graphics;
GA.Sim.EndCond = 2; % Run until converge (or fall)

% Set up the compass biped model
% GA.Sim.Mod = GA.Sim.Mod.Set('damp',0.3);
GA.Sim.Mod = GA.Sim.Mod.Set('I',0,'damp',0,'A2T',0.16,'A2H',0.12);

% Set up the terrain
start_slope = 0;
GA.Sim.Env = GA.Sim.Env.Set('Type','inc','start_slope',start_slope);

% Initialize the controller
switch whichCPG
        case {'ConSpitz','ConSpitz2','ConSpitz3','ConSpitz_eq'}
        GA.Sim.Con = ConSpitz;
        GA.Sim.Con.MinSat = [-maxAnkle,-maxHip];
        GA.Sim.Con.MaxSat = [ maxAnkle, maxHip];
        GA.Sim.Con.stDim = 1;
        GA.Sim.Con.FAM = whichCPG;  
    case {'2_level_CPG','2_level_CPG2','2_level_CPG3'}
        GA.Sim.Con = Controller;
        GA.Sim.Con.MinSat = [-maxAnkle,-maxHip];
        GA.Sim.Con.MaxSat = [ maxAnkle, maxHip];
        GA.Sim.Con.nPairs = N;
        GA.Sim.Con.stDim = 4*N+1;
    otherwise
        GA.Sim.Con = Matsuoka;
        GA.Sim.Con.startup_t = 1.0; % Give some time for the neurons to converge
        % before applying a torque
        GA.Sim.Con.FBType = 0; % no slope feedback
        GA.Sim.Con.nPulses = N;
        GA.Sim.Con.stDim = 4*N;
        GA.Sim.Con = GA.Sim.Con.SetOutMatrix([nAnkle1,nAnkle1,nHip]);
        GA.Sim.Con.MinSat = [-maxAnkle,-maxHip];
        GA.Sim.Con.MaxSat = [ maxAnkle, maxHip];
end
% Simulation parameters
GA.Sim.IC = [start_slope, start_slope, 0, 0, zeros(1, GA.Sim.Con.stDim)];
dur = 20;
GA.Sim = GA.Sim.SetTime(0,0.03,dur);

% Some more simulation initialization
GA.Sim.Mod.LegShift = GA.Sim.Mod.Clearance;

% GA.FitFcn = {1, @MOOGA.VelFit;  % walking velocity
%              2, @MOOGA.NrgEffFit; % energy spent of a step (1/(1+5*COT))
%              3:10, @MOOGA.VelRangeFit; range of velocities the controller
%              manages to walk in - TIME CONSUMING, runs the simulation
%              many times at slowly rising velocities.
%              11, @MOOGA.EigenFit}; % Eigenvalues of discrete poincare map
% GA.FitIDs = [1,2,3]; % Velocity and average COT
% GA.FitMinMax = [1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1];

GA.FitFcn = {1, @MOOGA.ASVelFit; % new walking velocity function, measured from some (nominal - 5th) stance switch to last stance switch
             2, @MOOGA.NrgEffFit; % energy spent of a step (1/(1+5*COT))
             3, @MOOGA.EigenFit; % Eigenvalues of discrete poincare map
             4, @MOOGA.VelFit; % walking velocity
             5, @MOOGA.StoFit; % stochastic walking fitness
             6, @MOOGA.UpDownFit}; % min/max slopes
GA.FitIDs = [1,5,6]; % The fitness functions that are actualy being optimized, out of those calculated 
GA.FitMinMax = ones(1,6); % sign of optimizing (1 - maximize, -1 - minimized)

% setup for StoFit:
GA.TerVar = 1e-5;
GA.SF_xend = 3;
GA.nTerForSto = 10;
if exist(['sto_ter_array_fixed' num2str(runcount) '.mat']) == 2
    load(['sto_ter_array_fixed' num2str(runcount) '.mat'],'ppv','dppv');
    GA.ppv = ppv;
    GA.dppv = dppv;
else
[GA.ppv,GA.dppv] = Sto_Ter_Array_Gen(runcount,GA.nTerForSto,GA.TerVar,GA.SF_xend);
end

GA.NFit = size(GA.FitIDs,2);
GA.Sim.PMFull = 1; % Run poincare map on all coords

GA = GA.InitGen();

% Update MOOGA parameters after each generation
function GA = GenFcn(GA)
    % Increase allowed genome range
    % Interpolate range
%         i = min(GA.Progress,10)/10;
%         Range = i*Range2 + (1-i)*Range1;
%         GA.Gen.Range = Range;

    % Reduce mutation range
%         j = min(GA.Progress,30)/30;
    j = GA.Progress/GA.Generations;
    GA.Gen.MutDelta = (1-j)*MutDelta0 + MutDelta1*j;
end
GA.GenerationFcn = @GenFcn;

% TrySim = deepcopy(GA.Sim);
% TrySim = GA.Gen.Decode(TrySim,GA.Seqs(1,:,1));
% TrySim = TrySim.SetTime(0,0.01,10);
% TrySim = TrySim.Init();
% TrySim.Run()

GA = GA.Run();
GA.Plot('Fit');

end