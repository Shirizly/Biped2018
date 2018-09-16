function Sim = Sim_setup(seq)
    Dur = 'lc';%20;
    timestep = 0.01;
    

     load('MatsuokaGenome.mat','Keys','Range','N',...
         'nAnkle1','nHip','maxAnkle', 'maxHip','MutDelta0','MutDelta1');

    Gen = Genome(Keys, Range);

Sim = Simulation();
    Sim.Graphics = 1;
    Sim.EndCond = 2; % Run until converge (or fall)

    % Set up the compass biped model
    % GA.Sim.Mod = GA.Sim.Mod.Set('damp',0.3);
    Sim.Mod = Sim.Mod.Set('I',0,'damp',0,'A2T',0.16,'A2H',0.12);

    % Set up the terrain
    start_slope = 0;
    Sim.Env = Sim.Env.Set('Type','inc','start_slope',start_slope);

    % Initialize the controller
%     Sim.Con = Matsuoka;
%     Sim.Con.startup_t = 1.0; % Give some time for the neurons to converge
%     % before applying a torque
%     Sim.Con.FBType = 0; % no slope feedback
%     Sim.Con.nPulses = N;
%     Sim.Con.stDim = 4*N;
%     Sim.Con = Sim.Con.SetOutMatrix([nAnkle,nHip]);
%     Sim.Con = Sim.Con.SetAnkles_selection(2); % set the simulation to work with two ankle joints
%     Sim.Con.MinSat = [-maxAnkle,-maxHip];
%     Sim.Con.MaxSat = [ maxAnkle, maxHip];

Sim.Con = Controller;
Sim.Con.MinSat = [-maxAnkle,-maxHip];
Sim.Con.MaxSat = [ maxAnkle, maxHip];
Sim.Con.nPairs = N;
Sim.Con.stDim = 4*N+1;

        
    % Simulation parameters
    Sim.IC = [start_slope, start_slope, 0, 0, zeros(1, Sim.Con.stDim)];

    if isnumeric(Dur)
        tend = Dur;
        Sim.EndCond = 0;
        Sim.doGoNoGo = 0;
    else
        switch lower(Dur)
            case 'step'
                Sim.EndCond = [1,1];
                tend = 20;
            case {'lc','limit cycle','limitcycle','limit_cycle'}
                Sim.EndCond = 0;
                tend = 'inf';
%                 Sim.minDiff = 1e-5;
        end
    end

    Sim.Graphics = 1;



    % Simulation parameters
    Sim = Sim.SetTime(0,timestep,tend);


    % decode the controller
    Sim = Gen.Decode(Sim, seq);
    % Set internal parameters (state dimensions, events, etc)
    Sim = Sim.Init();
    
    % Some more simulation initialization
    Sim.Mod.LegShift = Sim.Mod.Clearance;
    Sim.Con = Sim.Con.HandleEvent(1, Sim.IC(Sim.ConCo));
%     Sim.Con = Sim.Con.Adaptation();
end