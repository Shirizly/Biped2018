function [ sim ] = Run( sim )
% Run the simulation until an event occurs
% Handle the event and keep running
    X = [];
    sim.Out.T = [];
    Torques = [];
    Slopes = [];
    
    last_t = 1;
    
    if sim.Graphics == 1
        options=odeset('MaxStep',sim.tstep/10,'RelTol',1e-8,'AbsTol',1e-9,...
            'OutputFcn', @sim.Render, 'Events', @sim.Events);
    else
        options=odeset('MaxStep',sim.tstep/10,'RelTol',1e-8,'AbsTol',1e-9,...
            'Events', @sim.Events);
    end
    
    if sim.tstep ~= 0.0111
        tspan = sim.tstart:sim.tstep:sim.tend;
    else
        tspan = [sim.tstart,sim.tend];
    end
    [TTemp,XTemp,TE,YE,IE] = ...
        ode45(@sim.Derivative,tspan,sim.IC,options); %#ok<ASGLU>
    
    if sim.infTime == 1
        TimeCond = true;
        sim.tend = sim.tend + TTemp(end)-tspan(1);
    else
        TimeCond = TTemp(end)<sim.tend;
    end
    
    % Save state and time
    SuppPos = [ones(length(TTemp),1)*sim.Mod.xS, ...
               ones(length(TTemp),1)*sim.Mod.yS];
    X = [X; XTemp];
    sim.Out.T = [sim.Out.T; TTemp];
    
    if sim.nOuts>0
        % Save torques & slope
        Torques = sim.Con.Output(TTemp, ...
            XTemp(:,sim.ConCo)', XTemp(:,sim.ModCo)')';
        ThisSlope = sim.Env.SurfSlope(sim.Mod.xS);
        Slopes = repmat(ThisSlope,length(TTemp),1);
    end

    while TimeCond && sim.StopSim == 0
        % Deal with it
        StoreIC = 0;
        Xa = XTemp(end,:);
        for ev = 1:length(IE)
            % Is it a model event?
            ModEvID = find(IE(ev) == sim.ModEv,1,'first');
            if ~isempty(ModEvID)
                if ModEvID == 1 % Ground contact
                    % Check contact
                    [xNS, yNS]=sim.Mod.GetPos(Xa(sim.ModCo),'NS');
                    if abs(yNS-sim.Env.Surf(xNS))>1e-3
                        continue
                    end
                end
                        
                [sim.Mod,Xa(sim.ModCo)] = ...
                    sim.Mod.HandleEvent(ModEvID, ...
                        XTemp(end,sim.ModCo),TTemp(end));

                % Handle event interactions
                if ModEvID == 1 % Ground contact
                    StoreIC = 1; % Store the initial conditions right after impact

                    [sim.Con, Xa(sim.ModCo), Xa(sim.ConCo)] = ...
                        sim.Con.HandleExtFB(Xa(sim.ModCo),...
                        Xa(sim.ConCo),sim.Env.SurfSlope(sim.Mod.xS));

                    sim = sim.UpdateStats(TTemp,XTemp);
                    
                    % Shorten new swing leg (for ground clearance)
                    sim.Mod.LegShift = sim.Mod.Clearance;

                    if ~ischar(sim.Mod.curSpeed)
                        sim.TimeStr = ['t = %.2f s\nOsc.=%.3f\n',...
                         'Slope = %.2f ',char(176)','\nSpeed = %.3f m/s'];
                    end
                end
                
                if ModEvID == 2 % Robot fell down (hip height too low)
                    sim.Out.Type = 1;
                    sim.Out.Text = 'Robot fell down (hip height too low)';
                    sim.StopSim = 1;
                    break;
                end
            end

            % Is it a controller event?
            ConEvID = find(IE(ev) == sim.ConEv,1,'first');
            if ~isempty(ConEvID)
                [sim.Con,Xa(sim.ConCo)] = ...
                    sim.Con.HandleEvent(ConEvID, XTemp(end,sim.ConCo));

                % Handle event interactions
                switch ConEvID
                    case 1 % Neuron fired
%                         sim.Mod.LegShift = sim.Mod.Clearance;
                    case 2 % Leg extension
                        sim.Mod.LegShift = 0;
                end 
            end
        end

        if sim.nOuts>0 && sim.EndZMP == 1
            % Check ZMP
            ThisTorques = sim.Con.Output(TTemp, ...
                XTemp(:,sim.ConCo)', XTemp(:,sim.ModCo)')';
            CleanTorques = sim.GetCleanPulses(TTemp,XTemp,ThisTorques);
            Slope = sim.Env.SurfSlope(sim.Mod.xS);
            [ZMPfront, ZMPback] = sim.Mod.GetZMP(XTemp,CleanTorques,Slope);
            if ZMPfront>sim.Mod.A2T || abs(ZMPback)>sim.Mod.A2H
                % ZMP crossed the limits
                sim.Out.Type = 8;
                sim.Out.Text = 'ZMP crossed the limits';
                sim.StopSim = 1;
                break;
            end                
        end
        
        % Check ground clearance
%         [xNS,yNS] = sim.Mod.GetPos(Xa(sim.ModCo),'NS');
%         if yNS-sim.Env.Surf(xNS)<-1e-4*sim.Mod.L
%             if sim.Mod.LegShift>0
%                 % Robot hit the ground before extending the leg
%                 sim.Out.Type = 3;
%                 sim.Out.Text = 'Robot hit the ground before extending the leg';
%                 sim.StopSim = 1;
%             else
%                 if ~any(IE == 1)
%                     % Call impact handlers
%                     [sim.Mod,Xa(sim.ModCo)] = ...
%                         sim.Mod.HandleEvent(1, Xa(sim.ModCo),TTemp(end));
% 
%                     % Handle event interactions
%                     StoreIC = 1; % Store the initial conditions right after impact
% 
%                     [sim.Con, Xa(sim.ModCo), Xa(sim.ConCo)] = ...
%                         sim.Con.HandleExtFB(Xa(sim.ModCo),...
%                             Xa(sim.ConCo),sim.Env.SurfSlope(sim.Mod.xS));
% 
%                     sim = sim.UpdateStats(TTemp,XTemp);
% 
%                     if ~ischar(sim.Mod.curSpeed)
%                         sim.TimeStr = ['t = %.2f s\nOsc.=%.3f\n',...
%                          'Slope = %.2f ',char(176)','\nSpeed = %.3f m/s'];
%                     end
%                 else
%                     % Robot fell down
%                     sim.Out.Type = 1;
%                     sim.Out.Text = 'Robot fell down (leg pierced floor)';
%                     sim.StopSim = 1;
%                     break;
%                 end
%             end
%         end
        
        % Set new initial conditions
        sim.IC = Xa;
        if StoreIC
%             disp(['Time: ',num2str(TE)]);
            sim.ICstore(:,2:end) = sim.ICstore(:,1:end-1);
            sim.ICstore(:,1) = sim.IC';
            sim = sim.CheckConvergence();
        end
        
        if sim.StopSim
            break;
        end

        % Check for runaway events
        if any(abs(sim.IC(sim.ModCo(1:2)))>2*pi/3)
            % Leg angles are above 120 degrees 
            sim.StopSim = 1;
            break;
        end
        
        % Continue simulation
        if sim.tstep ~= 0.0111
            tspan = TTemp(end):sim.tstep:sim.tend;
        else
            tspan = [TTemp(end),sim.tend];
        end
        if length(tspan)<2
            % Can happen at the end of tspan
            break;
        end
        [TTemp,XTemp,TE,YE,IE] = ...
            ode45(@sim.Derivative,tspan,sim.IC,options); %#ok<ASGLU>
        
        if sim.infTime == 1
            TimeCond = true;
            sim.tend = sim.tend + TTemp(end)-tspan(1);
        else
            TimeCond = TTemp(end)<sim.tend;
        end
        
        % Save state and time
        SuppPos = [SuppPos;
                   ones(length(TTemp),1)*sim.Mod.xS, ...
                   ones(length(TTemp),1)*sim.Mod.yS]; %#ok<AGROW>
        X = [X; XTemp]; %#ok<AGROW>
        sim.Out.T = [sim.Out.T; TTemp];
        
        if sim.nOuts>0
            % Save torques & slope
            ThisTorques = sim.Con.Output(TTemp, ...
                XTemp(:,sim.ConCo)', XTemp(:,sim.ModCo)')';
            ThisSlope = sim.Env.SurfSlope(sim.Mod.xS);
            Torques = [Torques; %#ok<AGROW>
                       ThisTorques];
            Slopes = [Slopes; %#ok<AGROW>
                      repmat(ThisSlope,length(TTemp),1)];
                  
            if StoreIC && sim.CheckTorque == 1
                new_t = length(sim.Out.T);
                % Check that torque isn't 0 or constant
                if all(mean(abs(Torques(last_t:new_t,:)),1) < 1e-2)
                    % Torque output is 0 for the whole step
                    sim.Out.Type = 9;
                    sim.Out.Text = 'Torque signal is zero';
                    sim.StopSim = 1;
                    break
                end
                if all(std(Torques(last_t:new_t,:),1) < 1e-2)
                    % Torque output is constant for the whole step
                    sim.Out.Type = 9;
                    sim.Out.Text = 'Torque signal is constant';
                    sim.StopSim = 1;
                    break
                end
                
%                 figure
%                 plot(sim.Out.T(last_t:new_t), ...
%                      Torques(last_t:new_t,:))
%                 figure(1)
%                 disp(mean(abs(Torques(last_t:new_t,:)),1))
%                 disp(std(Torques(last_t:new_t,:),1))
%                 disp([sim.Con.tau, sim.Con.Amp']);

                last_t = new_t+1;
            end
        end
    end
    
    % Prepare simulation output
    sim.Out.X = X;
    if ~isempty(sim.Period)
        sim.Out.Tend = sim.Out.T(end);
    else
        sim.Out.Tend = sim.tend;
    end
    sim.Out.SuppPos = SuppPos;
    sim.Out.Torques = Torques;
    sim.Out.Slopes = Slopes;
    sim.Out.nSteps = sim.StepsTaken;
    sim.Out.StepsSS = sim.stepsSS;
    sim.Out.MaxSlope = [sim.MinSlope, sim.MaxSlope];
    
    if sim.Out.Type == 5
        % Simulation converged, calculate walking period
        sim.Period = sim.GetPeriod(1);
    end
end

