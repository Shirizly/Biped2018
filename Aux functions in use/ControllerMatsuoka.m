classdef ControllerMatsuoka < handle & matlab.mixin.Copyable
    
    properties
        %         CycleTime = [];
        %         Amp = [];
        %         Phase = [];
        
        X0 = []; % Membrane potential vector
        V0 = []; % Self-inhibition vector
        tau = []; % Fast time scaling parameter
        T = []; % Slow Time scaling parameter
        a = []; % Sinaptic strength parameter
        b = []; % Adaptation parameter
        tin = []; % Tonic input
        
        A = []; % Sinaptic strength matrix
        B = []; % Adaptation matrix
        Tin = []; % Tonic Input matrix
        
        IC;
        
        
    end
    
    methods
        % %%%%%% % Class constructor % %%%%%% %
        function C = ControllerMatsuoka(varargin)
            
            switch nargin
                case 0
                    C.X0 = zeros(2,1);
                    C.V0 = zeros(2,1);
                    C.tau = 1;
                    C.T = 100;
                    C.a = 1.5;
                    C.b = 1;
                    
                case 2
                    C.a = varargin{1};
                    C.b = varargin{2};
                case 4
                    C.tau = varargin{1};
                    C.T = C.tau*12;
                    C.a = varargin{4};
                    C.b = varargin{2};
                    C.tin = varargin{3};
                    C.IC = [0.33*C.tin*ones(1,2)* ...
                        diag(repmat([1, 0],1,1)), ...
                        zeros(1, 2)];
                    
                    C.X0 = C.IC(1:2).';
                    C.V0 = zeros(2,1);
                    
            end
            C.A = -(eye(2)-ones(2))*C.a;
            C.B = eye(2)*C.b;
            C.Tin = ones(2,1)*C.tin;
        end
        
        function Torques = Output(C, t) %
            
            Tspan = 0:0.001:t;
            Q0 = [C.X0;C.V0];
            
            [~,Q] = ode45(@C.Derivative,Tspan,Q0);
            
            X = Q(1:2,:);
            V = Q(3:4,:);
            
            Y = poslin(X);
            Torques = Y(1)-Y(2);
            
            
            
            
        end
        function [dQ] = Derivative(C, t, Q) % Defining the equations of the neurons to be integrated
            X = Q(1:2);
            V = Q(3:4);
            
            Y = max(X,0);
            
            dX = (C.Tin-C.A*Y-C.B*V-X)/C.tau;
            dV = (Y-V)/C.T;
            dQ = [dX;dV];
        end
    end
end