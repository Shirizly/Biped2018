function [T,R] = STRcheck2(varargin)
%Stochastic Terrain Robustness Check:
%The function checks how the model and controller are affected by a
%stochastic ground. Each robot needs to walk N steps on the stochastic
%terrain to confirm his robustness. For each amplitude of terrain, S robots
%are checked for robustness. 
%The function returns vector T of Terrain amplitude, and vector R of
%percentage of robust robots in the correlating Terrain amplitude.

%Start variables:
%(Con,CB,x_nom,N,S)

%INPUTS:
%Con- Controller
%CB- Model
%x_nom- fixed point for limit cycle
%N - Number of steps to confirm robustness
%S - Number of robots to check for each Terrain amplitude.
Graphics = 0; %Default

switch nargin
    case 1 %open backuplog to resume calculations
        load(varargin{1});
        Graphics = 0;
    case 2 %open backuplog to resume calculations+Graphics
        load(varargin{1});
        Graphics=varargin{2};
        display(varargin);
    case 5 %new calculation
        T = [0];
        R = [];
        dT = 0.005;
        Ttemp = T(1);
        stopflag=1;
        Con=varargin{1};
        CB=varargin{2};
        x_nom=varargin{3};
        N=varargin{4};
        S=varargin{5};
     case 6 %new calculation+Graphics
        T = [0];
        R = [];
        dT = 0.005;
        Ttemp = T(1);
        stopflag=1;
        Con=varargin{1};
        CB=varargin{2};
        x_nom=varargin{3};
        N=varargin{4};
        S=varargin{5};
        Graphics=varargin{6};
     case 7 %new calculation+Graphics
        T = [0];
        R = [];
       
        Ttemp = T(1);
        stopflag=1;
        Con=varargin{1};
        CB=varargin{2};
        x_nom=varargin{3};
        N=varargin{4};
        S=varargin{5};
        Graphics=varargin{6};
        dT=varargin{7};
end
clear varargin

    while(stopflag)
        display(['Checking Chaos Parameter ', num2str(Ttemp)]);
        cnt=0;
        for i=1:S
               [fun,fundot]=StochTerrain(Ttemp,N);
               Env = Terrain(5,fun,fundot); %Need to create a stochastic Terrain
               CB_temp=copy(CB);
               Sim = Simulation3(CB_temp,Con,Env);
               Sim.Graphics = 0;
                Sim.EndCond = 2;
                Sim.IC = x_nom;
                Sim = Sim.SetMaxSteps(N);

                %Some more simulation initialization
                Sim.Mod.LegShift = Sim.Mod.Clearance; 
                %%
                Sim = Sim.SetTime(0,0.01,3*N);
                % Set internal parameters (state dimensions, events, etc)
                %%
                Sim = Sim.Init();
                %%
                % Simulate
                Sim = Sim.Run();

                if(Sim.Out.Type==4)
                    cnt=cnt+1;
                    Approval='Success';
                else
                    Approval='Failure';
                end
                display(['Robot No. ' int2str(i),': ',Approval]);
                %close all;

        end
        R = [R , cnt];

        
        if(Graphics)
            if(length(T)>1)
                plot([T(end-1),T(end)],[(R(end-1)/S)*100,(R(end)/S)*100],'-or');
                shg;
            else
                figure(1);
                hold on;
                plot((T(end)/S)*100,R(end),'or');
                title('Stochastic Ground Robustness');
                xlabel('Variance');
                ylabel('Successful Percentage');
                drawnow;
                pause(1);
                refreshdata;
            end
            
        end
        
        if(cnt<=0.05*S)
            stopflag=0;
        else
            T=[ T , Ttemp+dT];
            Ttemp=Ttemp+dT;
        end
        
        save('BackupRobustnessLog2');
        pause(1);


    end









end