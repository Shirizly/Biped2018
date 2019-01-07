classdef ConSpitz < handle & matlab.mixin.Copyable

    
    properties
        name = 'Spitz CPG'
        
       
        FAM = 0; % frequency adaptation method
        % Slope feedback parametrs
		FBType = 0;  % 0 - no slope feedback
        
        % Saturation
        MinSat; MaxSat;
        
       
        % parameters for MOGA use:
        des_period = 1.4; % desired period for rescaling method.
        
        
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
        omega = 1.106512566;
        omega0 = 1.106512566;
        
        % Initial phase for oscillator:
        P_0 = 0;
                
        stDim = 1; % state dimension
        nEvents = []; % automatically updated in addPulse
        
        % Controller Output
        nPulses = 0; % Overall number of pulses - automatically updated in addPulse
        OutM;        % Output matrix (multiplies Switch vector)
        
        Amp = 0;     % Pulse amplitude in N
        Amp0 = 0;    % Base pulse amplitude in N
        Offset = 0;  % Defines beginning of pulse as % of neuron period
        Duration = 0;% Pulse duration as % of neuron period
        
        Switch;      % 0 when off, Amp when pulse is active
        startup_t;
                ExtPulses = []; % External pulses IDs
        % Feedback (on the PG level)
%         FBType = 2;  % 0 - no feedback
                     % 1 - single gain for omega and each joint
                     % 2 - individual gains for each pulse (not
                     % implemented?)
        lastPhi = 0;
%         s_in = 0;    % Speed input - higher level input to control the
                     % desired walking speed
         IC_MO = [0,1.106512566];
         startup_c;
         omega_c;
        % Phase reset
         ExtP_reset = 0; % set to a certain phase to use phase reset
                         % on foot contact - feedback reset
                         
        % Angular velocity impulses
        FBImpulse = 0; % 0 - no feedback
                       % 1 - add a constant value
                       % 2 - set ang. vel. to certain value
        AngVelImp = [];
                
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
        
%       Adaptation parameters:
        k_a = zeros(2,2); % matrix tying the feedback (angle sum and difference) to the parameters (tau and amp)
        dtnom = 0.3;
        dtconv = 1;
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

            % Apply saturation
            if ~isempty(NC.MinSat)
                for j = 1:length(NC.MinSat)
                    Torques(j,:) = ...
                        min(max(Torques(j,:),NC.MinSat(j)),NC.MaxSat(j));
                end
            end

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
        
         function [NC] = Adaptation(NC,X)
			% % % slope adaptaion:
            switch NC.FBType
                case 0
                    % NO FEEDBACK
                case 1
                    % Apply higher-level adaptation:
                    P0 = [NC.omega0;0];
                    FB = (X(1)+X(2))/2; % FB is: slope's angle
                    P = P0+NC.k_a*FB;
                    NC.omega = max(P(1),0.02);
                    NC.Amp = NC.Amp0+P(2)*[-1,1,1];
                case 2 % learning the nominal delta-theta for the adaptation
                    dt = X(1)-X(2);
                    NC.dtconv = abs(NC.dtnom - dt);
                    NC.dtnom = dt;
                case 3
                    % Apply higher-level adaptation:
                    P0 = [NC.omega0;0];
                    t1 = X(1);
                    t2 = X(2);
                    FB = [(t1+t2)/2;(t1-t2)-NC.dtnom]; % FB is: slope's angle, difference in relative angle between the legs (compared to nominal)
                    P = P0+NC.k_a*FB;
                    NC.omega = max(P(1),0.02);
                    NC.Amp = NC.Amp0+P(2)*[-1,1,1];
            end
        end
        
        % %%%%%% % Derivative % %%%%%% %
        function [Xdot] = Derivative(NC, ~, Xmod, X) % here MO stands for the controller
            Xdot = NC.omega;
        end
        
        % %%%%%% % Events % %%%%%% %
        function [value, isterminal, direction] = Events(NC, X)
            value = ones(NC.nEvents,1);
            isterminal = ones(NC.nEvents,1);
            direction = -ones(NC.nEvents,1);
            
%             Check for end of period
             value(1) = NC.P_th - X(end);

            % Check for leg extension signal (for clearance)
%             value(2) = NC.P_LegE - X;

            % Check for Neuron firing (E of hip)

            
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

                    
%                 case 2
%                     % reset the phase according to the hip E neuron
%                     for i=1:NC.nPulses
%                         NC.Switch(i) = 0; % switching off all pulses
%                         if NC.Offset(i)==0
%                             % switch on only pulses that start at the
%                             % begining of the step period
%                             NC.Switch(i) = NC.Amp(i);
%                         end
%                     end
%                     
%                     if NC.omega_c==0 %%learning the MO natural period, to match the PG's frequency
%                         NC.omega_c = 1;
%                     else
%                         if NC.startup_c == 0 || contains(NC.FAM,'2_level_CPG2')
%                         NC.omega = NC.omega/(Xa(end)); % this allows the PG to update on the MO frequency, either once or once per cycle
%                         NC.startup_c = 1;
% %                         disp('start pulses')
%                         else
%                             if contains(NC.FAM,'2_level_CPG3')
% %                             Xa(end) % Just here to see the phase difference in real time simulation
%                                 NC.omega = NC.omega*(1+ NC.k_om*(1-Xa(end))); % 3rd method frequency adaptation
%                             end
%                         end
%                         
%                         % need to add condition for doing this only the
%                         % first time,
%                         if NC.IC_MO(6) == 0
%                             NC.IC_MO = [Xb,NC.omega];
%                             NC.IC_MO(end-1) = 0;
%                         end
%                         
%                     end
%                     Xa(end) = NC.P_reset; % reset phase
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
                
                
                NC = NC.Adaptation(Xmod);
                
                for i=1:NC.nPulses
                    NC.Switch(i) = 0; % switching off all pulses
                    if NC.Offset(i)==0;
                        % Turn on signal now
                        NC.Switch(i) = NC.Amp(i);
                    end
                end
                Xcon = NC.ExtP_reset; % reset phase
            end

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