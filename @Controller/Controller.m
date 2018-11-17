classdef Controller < handle & matlab.mixin.Copyable

    
    properties
        name = '2 level CPG'
        
        %% Matsuoka parameters
        %         Parameters for the CPG
        tau0 = 1; tau = 1;
        tav0 = 1; tav = 1;
        tau_ratio = 12; 
        beta = 0.1;
        u0 = 1;
        win = [];
        wex = [];
        W = [];
     
        stDim = 4; % state dimension
%         nEvents = 1; % num. of simulation events
        
		% CPG feedback parameters:
		k_hip_fb = 0;
        legSwitch = 1; % switch the feedbcak sign when the legs are switching
		
		% Slope feedback parametrs
		FBType = 0;  % 0 - no slope feedback
		
        % Controller Output
        startup_t = 10; % torque start time. give the CPG time to converge 
        startup_c = 0; % torque start check, based on period start
        omega_c = 0;
        nPairs = 1; % number of flexor/extensor pairs which actuate the joint
        njoints = 2; % number of actuated joints in the CB
        OutM = [];
        NAmp0 = 1.5;
        NAmp = 1.5;
        ExtPulses = [];
        
        % Saturation
        MinSat; MaxSat;
        
        % Higher level speed command parameters
        s_in = 0;    % Speed input - higher level input to control the
                     % desired walking speed
        ks_tau = 0;	% gain on tau
        ks_out = 0; % gain on the tonics inputs
        
        % parameters for MOGA use:
        des_period = 1.4; % desired period for rescaling method.
        
        % Phase reset
        ExtP_reset = []; % set to a certain phase to use phase reset
                         % on foot contact
        
        % Set keys
        SetKeys = {'npulses', 'nneurons', 'n_pulses', 'n_neurons', ...
            'tau', 'tau_u', 'tav', 'tau_v', 'tau_r', '\tau_r', ...
            'tau_ratio', '\tau_ratio', ...
            'beta', 'win', 'wex', 'weights', 'amp0', 'amp', ...
            'ks_tau', 'speed_tau', 'tau_speed_gain', 'ks_\tau', ...
            'ks_out', 'speed_out', 'torque_speed_gain', 'ks_c',...
            '2neuron_symm_weights','2neuron_general_weights',...
            'amp_2n_dif_inputs','amp_2n_same_inputs','ks_c_2n_symm','ks_c_2n_general',...
            'amp_6n_symm','6neuron_taga_like',...
            '6neuron_taga_like_symm','k_hip_fb','des_period','k_om'};
        %% End Matsuoka parameters, start of 2nd level CPG parameters
        
        % LIF parameters
        P_reset = 0;
        P_th = 1;
        
        P_LegE = 0.65; % timed leg extension

        % Oscillator's frequency
        omega = 1; % 1.106512566;
        omega0 = 1; % 1.106512566;
        
        % Initial phase for oscillator:
        P_0 = 0;
%         startup_t = 0;  % Time before any torque output is available
                
%         stDim = 1; % state dimension
        nEvents = []; % automatically updated in addPulse
        
        % Controller Output
        nPulses = 0; % Overall number of pulses - automatically updated in addPulse
%         OutM;        % Output matrix (multiplies Switch vector)
        
        Amp = 0;     % Pulse amplitude in N
        Amp0 = 0;    % Base pulse amplitude in N
        Offset = 0;  % Defines beginning of pulse as % of neuron period
        Duration = 0;% Pulse duration as % of neuron period
        
        Switch;      % 0 when off, Amp when pulse is active
%         ExtPulses;   % External pulses IDs
        
        % Feedback (on the PG level)
%         FBType = 2;  % 0 - no feedback
                     % 1 - single gain for omega and each joint
                     % 2 - individual gains for each pulse (not
                     % implemented?)
        lastPhi = 0;
%         s_in = 0;    % Speed input - higher level input to control the
                     % desired walking speed
                     
        % Phase reset
%         ExtP_reset = []; % set to a certain phase to use phase reset
                         % on foot contact
                         
        % Angular velocity impulses
        FBImpulse = 0; % 0 - no feedback
                       % 1 - add a constant value
                       % 2 - set ang. vel. to certain value
        AngVelImp = [];
        
        IC_MO = zeros(1,6); % save and load IC for beginning of controller's period
        
        % Gains
        kOmega_u = 0.0;
        kOmega_d = 0.0;
        kTorques_u = 0;
        kTorques_d = 0;
        sOmega_f = 0.0;
        sOmega_s = 0.0;
        sTorques_f = 0;
        sTorques_s = 0;
        k_om = 0;
        
        % Saturation
%         MinSat; MaxSat;
        
        % Set keys
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function [NC] = Controller(varargin)
            % Set torque parameters
%             NC.nPulses = 4;
%             NC.Amp0=[-13.6255 -4.6281 13.6255 0];
            NC.nPulses = 0;
            NC.Amp0=[];
            NC.Amp=NC.Amp0;
%             NC.OutM = [1 1 0 0 ; 0 0 1 1]; % updated automatically by the
%             AddPulse function
            


            
            NC.Switch=zeros(NC.nPulses,1);
            
            % Set adaptation parameters
            NC.kOmega_u = 0; %0.4579;
            NC.kOmega_d = 0; %0.9;
            NC.kTorques_u= 0; %[95 -443 95 0];
            NC.kTorques_d= 0; %[80 -295 80 0];
            
            NC.nEvents=2+NC.nPulses*2;
        end
        
        function [NC] = ClearTorques(NC)
            NC.nPulses = 0;
            NC.OutM = [0,0]';
            NC.Amp0 = [];
            NC.Amp = [];
            NC.Offset = [];
            NC.Duration = [];
            NC.Switch = 0;
            if NC.FBType == 2
                NC.kTorques_u = [];
                NC.kTorques_d = [];
            end
            NC.nEvents = 1; 
            NC.ExtPulses = [];
        end
        
        function NC = Reset(NC,Phase) % checks which pulses should be active, and activate them
            if nargin<2
                NC.Switch = 0*NC.Switch;
            else
                % Find which Torques should be active
                Start = NC.Offset;
                if ~isempty(Start)
                    End = Start + NC.Duration;
                    [Xperc] = NC.GetPhasePerc(Phase);
                    On = Xperc >= Start & Xperc < End;
                    NC.Switch = (NC.Amp.*On)';
                else
                    NC.Switch = 0;
                end
            end 
        end
        
        function [Torques] = Output(NC, t, ~, ~)

            Torques = repmat(NC.OutM*NC.Switch,1,length(t));
%             Torques(:,t<NC.startup_t) = 0*Torques(:,t<NC.startup_t);
            % Apply saturation
            if ~isempty(NC.MinSat)
                for j = 1:length(NC.MinSat)
                    Torques(j,:) = ...
                        min(max(Torques(j,:),NC.MinSat(j)),NC.MaxSat(j));
                end
            end
            Torques = NC.startup_c*Torques;
           


        end

        function [per] = GetPeriod(NC)
            per = (NC.P_th-NC.P_reset)/NC.omega;
        end
        
        function [Xperc] = GetPhasePerc(NC,X)
            Xperc = (X-NC.P_reset)/(NC.P_th-NC.P_reset);
        end
        
        function diff = PhaseDiff(NC,ph1,ph2)
            diff = ph1 - ph2;
            % wrap it up
            Range = NC.P_th - NC.P_reset;
            diff(diff>Range/2) = diff(diff>Range/2) - Range;
            diff(diff<-Range/2) = diff(diff<-Range/2) + Range;
        end
        
        
        
        % %%%%%% % Derivative % %%%%%% %
        function [Xdot] = Derivative(NC, ~, Xmod, X) 
            phidot = NC.omega;
            
            % odd - flexor
            % even - extensor
            
            % % feed back to hip:
            feedE1 = zeros(2*NC.nPairs,1);
            feedE1(end-1:end,1) = [-NC.k_hip_fb * NC.legSwitch *(Xmod(2) - Xmod(1));...
                +NC.k_hip_fb * NC.legSwitch * (Xmod(2) - Xmod(1));];
%             feedE1(end-2:end-1,1) = [-NC.k_hip_fb * NC.legSwitch *(Xmod(2) - Xmod(1));... % this is the wrong, all the results up to 28.03 use this
%                 +NC.k_hip_fb * NC.legSwitch * (Xmod(2) - Xmod(1));];
            
            % % Maybe add synch to the hip joint angular speed?
			% feedE2 = [0;0;0;0;...
            %     -MO.k_d_hip_fb * MO.legSwitch *(Xmod(3) - Xmod(4));...
            %     +MO.k_d_hip_fb * MO.legSwitch * (Xmod(3) - Xmod(4));];
			
            % X = [u_i; v_i];
            u = X(1:2*NC.nPairs,1);
            v = X(2*NC.nPairs+1:end-1,1);
            y = max(u,0);
            udot = 1/NC.tau*(-u - NC.beta*v + NC.NAmp - NC.wex*y +feedE1);
            vdot = 1/NC.tav*(-v+y);
            
            Xdot = [udot;
                    vdot
                    phidot];
        end
        
        % %%%%%% % Events % %%%%%% %
        function [value, isterminal, direction] = Events(NC, X)
            value = ones(NC.nEvents,1);
            isterminal = ones(NC.nEvents,1);
            direction = -ones(NC.nEvents,1);
            
            % Check for end of period
%             value(1) = NC.P_th - X(end);

            % Check for leg extension signal (for clearance)
%             value(2) = NC.P_LegE - X;

            % Check for Neuron firing (E of hip)
            u = X(1:2*NC.nPairs,1);
            value(2) = u(end-1)-u(end);
            direction(2) = 1;
            
            % Check for switching on/off signal
            Xperc = ones(NC.nPulses,1)*NC.GetPhasePerc(X(end));
            value(3:2+NC.nPulses,1) = NC.Offset.' - Xperc;
            value(3+NC.nPulses:2+2*NC.nPulses,1) = ...
                NC.Offset.' + NC.Duration.' - Xperc;
        end
        
        function [NC,Xa] = HandleEvent(NC, EvID, Xb)
            Xa = Xb;
            switch EvID
%                 case 1
%                     % original phase reset (phase equal to 1)
%                     for i=1:NC.nPulses
%                         NC.Switch(i) = 0; % switching off all pulses
%                         if NC.Offset(i)==0;
%                             % Turn on signal now
%                             NC.Switch(i) = NC.Amp(i);
%                         end
%                     end
%                     Xa(end) = NC.P_reset; % reset phase

                    
                case 2
                    % reset the phase according to the hip E neuron
                    for i=1:NC.nPulses
                        NC.Switch(i) = 0; % switching off all pulses
                        if NC.Offset(i)==0
                            % switch on only pulses that start at the
                            % begining of the step period
                            NC.Switch(i) = NC.Amp(i);
                        end
                    end
                    
                    if NC.omega_c==0 %%learning the MO natural period, to match the PG's frequency
                        NC.omega_c = 1;
                    else
                        if NC.startup_c == 0 || NC.k_om == 0
                        NC.omega = NC.omega/(Xa(end)); % this allows the PG to update on the MO frequency, either once or once per cycle
                        NC.startup_c = 1;
%                         disp('start pulses')
                        else
%                             Xa(end) % Just here to see the phase difference in real time simulation
                            NC.omega = NC.omega*(1+ NC.k_om*(1-Xa(end))); % 3rd method frequency adaptation
                        end
                        
                        % need to add condition for doing this only the
                        % first time,
                        if NC.IC_MO(6) == 0
                            NC.IC_MO = [Xb,NC.omega];
                            NC.IC_MO(end-1) = 0;
                        end
                        
                    end
                    Xa(end) = NC.P_reset; % reset phase
%                     disp('phase reset')

                case num2cell(3:2+NC.nPulses)
                    % Switch on signal
                    NC.Switch(EvID-2) = NC.Amp(EvID-2);
                case num2cell(3+NC.nPulses:2+2*NC.nPulses)
                    % Switch off signal
                    PulseID = EvID-(2+NC.nPulses);
                    NC.Switch(PulseID) = 0;
                    if any(PulseID == NC.ExtPulses)
                        % Set offset back to 200%
                        NC.Offset(PulseID) = 2;
                    end
            end
        end
        
        function [NC, Xmod, Xcon] = HandleExtFB(NC, Xmod, Xcon, Slope)
            % This function is called when the leg hits the ground
            
            if ~isempty(NC.ExtP_reset)
                % Perform a phase reset
                Xcon = NC.ExtP_reset;
                
                % Check if any event happens at ExtP_reset
                [value, it, dir] = NC.Events(Xcon); %#ok<NASGU,ASGLU>
                EvIDs = find(value == 0);
                for ev = 1:length(EvIDs)
                    [NC,Xcon] = NC.HandleEvent(EvIDs(ev),Xcon);
                end
            end
            
            switch NC.FBImpulse
                case 1 % 1 - add a constant value
                    Xmod(3:4) = Xmod(3:4) + NC.AngVelImp;
                case 2 % 2 - set ang. vel. to certain value
%                     delta = NC.AngVelImp - Xmod(3:4)
                    Xmod(3:4) = NC.AngVelImp;
            end
            
            % Activate external pulses
%             NC.Switch(NC.ExtPulses) = NC.Amp(NC.ExtPulses);
%             NC.Offset(NC.ExtPulses) = NC.GetPhasePerc(Xcon);
%             % Set the torque to get turned off after the neuron fires if
%             % the off event is larger than 100% of osc. period
%             Overflow = NC.Offset(NC.ExtPulses)+NC.Duration(NC.ExtPulses)>1;
%             NC.Offset(NC.ExtPulses) = NC.Offset(NC.ExtPulses) - Overflow;
        end
        
        function PlotTorques(NC,tstep,mux)
            if nargin<2
                [Time,TorqueSig]=NC.GetTorqueSig();
            else
                if nargin<3
                    [Time,TorqueSig]=NC.GetTorqueSig(tstep);
                else
                    [Time,TorqueSig]=NC.GetTorqueSig(tstep,mux);
                end
            end

            plot(Time,TorqueSig);
        end
        
        function [Time,TorqueSig]=GetTorqueSig(NC,pstep,mux)
            % Check number of inputs to method
            if nargin<2
                pstep = 0.01;
            end
            if nargin<3
                mux = 1;
            end
            
            % Prepare empty variables
            Phase = (0:pstep:1)';
            Np = length(Phase);
            Time = Phase/NC.omega;
            if mux
                TorqueSig = zeros(Np,size(NC.OutM,1));
            else
                TorqueSig = zeros(Np,size(NC.OutM,2));
            end
            
            % Build torque signals
            for p = 1:NC.nPulses
                pStart = NC.Offset(p);
                pEnd = pStart + NC.Duration(p);
                if mux
                    j = find(NC.OutM(:,p)==1,1,'first');
                    TorqueSig(:,j) = TorqueSig(:,j) + ...
                        NC.Amp(p)*(Phase >= pStart & Phase < pEnd);
                else
                    TorqueSig(:,p) = ...
                        NC.Amp(p)*(Phase >= pStart & Phase < pEnd);
                end
            end
        end
        
        function diffs = CompareSigs(NC,NC2,pstep)
            if nargin<3
                pstep = 0.001;
            end
            [thisT,thisSig]=NC.GetTorqueSig(pstep,1);
            [~,otherSig]=NC2.GetTorqueSig(pstep,1);
            
            diffs = zeros(1,size(thisSig,2));
            try
                dSig = abs(thisSig-otherSig);
                for p = 1:size(thisSig,2)
                    diffs(p) = trapz(thisT,dSig(:,p));
                end
            catch err
                diffs = diffs + 1e3;
                disp(['Error comparing signals: ',err]);
            end
        end
            
    end
end