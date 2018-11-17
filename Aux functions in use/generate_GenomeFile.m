function generate_GenomeFile(whichCase)
%GENERATE_GENOMEFILE simple and easy function to generates the genome files.
% all in one place. so different sripts and function wouldn't overwrite the
% genomeFile.
% 
% Input:
% *) 'whichCase' - can be:
%       #) '2N_symm' - for 2neuron symmetric CPG
%       #) '2N_general' - for 2neuron general CPG
%       #) '4N_tagaLike' - for 4neuron taga-like CPG
%       #) '4N_general' - for 4neuron general CPG
% 
% NOTE: change the "keys" if you need to change the GA optimization
% parameters.

% genome file name:
genome_file = 'MatsuokaGenome.mat';

% define Mutation strength:
MutDelta0 = 0.2;   MutDelta1 = 0.05;

switch whichCase
    case {'2N_symm','2N_general'}
        nAnkle1 = 0; % Number of ankle1 torque pulses
        nAnkle2 = 0; % Number of ankle2 torques pulses
        nHip = 1;   % Number of hip torques
        % % Parameters:
        maxAnkle = 0;   % Max ankle torque
        maxHip = 20;    % Max hip torque
        maxW = 5;
        minW = 0;
    case {'6N_tagaLike_2Ank_torques','6N_tagaLike_2Ank_torques_symm',...
            '6N_tagaLike_2Ank_torques_symm_feedback',...
            '6N_tagaLike_2Ank_torques_symm_feedback_eq',...
            '6N_tagaLike_2Ank_torques_symm_feedback_eq_adaptation',...
            '6N_tagaLike_2Ank_torques_symm_with_rescale',...
            '6N_general_2Ank_torques','6N_general_2Ank_torques_feedback',...
            '2_level_CPG','2_level_CPG_eq','2_level_CPG2','2_level_CPG3',...
            'ConSpitz','ConSpitz2','ConSpitz3','ConSpitz_eq','ConSpitz_eq_adaptation'}
        nAnkle1 = 1; % Number of ankle1 torque pulses
        nAnkle2 = 1; % Number of ankle2 torques pulses
        nHip = 1;   % Number of hip torques
        % % Parameters:
        maxAnkle = 20;   % Max ankle torque
        maxHip = 20;    % Max hip torque
        maxW = 5;
        minW = 0;
    otherwise
        error('invalid input');
end

N = nAnkle1+nAnkle2+nHip;

disp('The defined Genome is:')
disp(['case: ',whichCase]);

switch whichCase
    case '2N_symm'
        Mw = maxW*[1,1];
        mw = 0*Mw;
        
        Keys = {'\tau_r', 'beta','amp_2n_same_inputs',    '2neuron_symm_weights', 'ks_\tau',     'ks_c_2n_symm', 'IC_matsuoka';
                      1 ,      1,                   2,                         1,        1 ,          1,            0 };
        Range = {  0.02 ,    0.2,                  mw,                         0,      -10 ,         -1; % Min
                   0.25 ,    2.5,                  Mw,                        10,       10 ,          1}; % Max

               disp('Tau_retio = 12');

        % Note: because of some old mistakes. the tonic input gene ('c') 
        %   is encoded with two values. only one of them is used.          
    case '2N_general'
        Mamp = maxHip*[1,1];
        mamp = 0*Mamp;
        
        Mw = maxW*[1,1];
        mw = 0*Mw;
        
        Keys = {'\tau_r', 'beta', 'amp_2n_dif_inputs',    '2neuron_general_weights', 'ks_\tau', 'ks_c_2n_general', 'IC_matsuoka';
                      1 ,      1,                   2,                            2,        1 ,                 2,            0 };
        Range = {  0.02 ,    0.2,                mamp,                           mw,      -10 ,           [-1,-1]; % Min
                   0.25 ,    2.5,                Mamp,                           Mw,       10 ,             [1,1]}; % Max
               disp('Tau_retio = 12');

    case '6N_tagaLike_2Ank_torques'
        
        Mamp = [maxAnkle*ones(1,2*nAnkle1),...
            maxAnkle*ones(1,2*nAnkle2),...
            maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        mw = 0*ones(1,6);
        Mw = maxW*ones(1,6);
        
        % CPG strucute: (ALSO Symm W_ij = W_ji)
        % 1) every pair of Extensor and reflexor neurons are connected to
        %   each other with symmetric weights.
        % 2) both of the ankles' Extensor neurons are connected to both of the hip neurons
        %
        %  A2_F    A2_E
        %(3) O-----O (4)
        %       /  |
        %      /   |
        %     /    |
        %    /     |
        %   H_F   H_E
        %(5) O-----O (6)  
        %     \    |             %   w = [0  , W12, 0  , 0  , 0  , 0; 
        %      \   |             %        W21, 0  , 0  , 0  , W25, W26;
        %       \  |             %        0  , 0  , 0  , W34, 0  , 0;
        %        \ |             %        0  , 0  , W43, 0  , W45, W46;
        %(1) O-----O (2)         %        0  , 0  , 0  , 0  , 0  , W56;
        %   A1_F    A1_E         %        0  , 0  , 0  , 0  , W65, 0;
        %                        
        %                        w12 = w21 = w34 = w43 = w1  
        %                        w56 = 65  = w2
        %                        w25 = w3
        %                        w26 = w4
        %                        w45 = w5
        %                        w46 = w6
        
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        Keys = {'\tau_r', 'beta', 'amp_6n_symm',   '6neuron_taga_like', 'ks_\tau',     'ks_c_6n_symm', 'IC_matsuoka';
                      1 ,      1,             1,                     6,        1 ,                 1 ,            0 };
        Range = {  0.02 ,    0.2,             0,                    mw,      -10 ,                 -0.1*maxAnkle; % Min
                   0.25 ,    2.5,      maxAnkle,                    Mw,       10 ,                  0.1*maxAnkle}; % Max
		disp('Tau_retio = 12');

	case '6N_tagaLike_2Ank_torques_symm' % Rea's main controller for his Thesis results

        Mamp = [maxAnkle*ones(1,2*nAnkle1),...
            maxAnkle*ones(1,2*nAnkle2),...
            maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        mw = 0*ones(1,4);
        Mw = maxW*ones(1,4);
        
        % CPG strucute: (ALSO Symm W_ij = W_ji)
        % 1) every pair of Extensor and reflexor neurons are connected to
        %   each other with symmetric weights.
        % 2) both of the ankles' Extensor neurons are connected to both of the hip neurons
        %
        %  A2_F    A2_E
        %(3) O-----O (4)
        %       /  |
        %      /   |
        %     /    |
        %    /     |
        %   H_F   H_E
        %(5) O-----O (6)  
        %     \    |             %   w = [0  , W12, 0  , 0  , 0  , 0; 
        %      \   |             %        W21, 0  , 0  , 0  , W25, W26;
        %       \  |             %        0  , 0  , 0  , W34, 0  , 0;
        %        \ |             %        0  , 0  , W43, 0  , W45, W46;
        %(1) O-----O (2)         %        0  , 0  , 0  , 0  , 0  , W56;
        %   A1_F    A1_E         %        0  , 0  , 0  , 0  , W65, 0;
        %                        
        %                        w12 = w21 = w34 = w43 = w1  
        %                        w56 = w65  = w2
        %                        w25 = w46 = w3
        %                        w26 = w45 = w4
        
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        Keys = {'\tau_r', 'beta', 'amp_6n_symm',  '6neuron_taga_like_symm', 'ks_\tau',     'ks_c_6n_symm', 'IC_matsuoka';
                      1 ,      1,             1,                         4,        1 ,                 1 ,            0 };
        Range = {  0.02 ,    0.2,             0,                        mw,      -10 ,                 -0.1*maxAnkle; % Min
                   0.25 ,    2.5,      maxAnkle,                        Mw,       10 ,                  0.1*maxAnkle}; % Max
               disp('Tau_retio = 12');
               
    case '6N_tagaLike_2Ank_torques_symm_with_rescale'
        Mamp = [maxAnkle*ones(1,2*nAnkle1),...
            maxAnkle*ones(1,2*nAnkle2),...
            maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        mw = 0*ones(1,4);
        Mw = maxW*ones(1,4);
        
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        Keys = {'\tau_r', 'beta', 'amp_6n_symm',  '6neuron_taga_like_symm', 'ks_\tau',     'ks_c_6n_symm', 'des_period', 'IC_matsuoka';
                      1 ,      1,             1,                         4,        1 ,                 1 ,            1,            0 };
        Range = {  0.02 ,    0.2,             0,                        mw,      -10 ,      -0.1*maxAnkle,            1; % Min
                   0.25 ,    2.5,      maxAnkle,                        Mw,       10 ,       0.1*maxAnkle,            2}; % Max
               disp('Tau_retio = 12');

	case '6N_tagaLike_2Ank_torques_symm_feedback'

        Mamp = [maxAnkle*ones(1,2*nAnkle1),...
            maxAnkle*ones(1,2*nAnkle2),...
            maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        mw = 0*ones(1,4);
        Mw = maxW*ones(1,4);
         
        Keys = {'\tau_r', 'beta', 'amp_6n_symm',  '6neuron_taga_like_symm', 'ks_\tau',     'ks_c_6n_symm','k_hip_fb', 'IC_matsuoka';
                      1 ,      1,             1,                         4,        1 ,                 1 ,         1,            0 };
        Range = {  0.02 ,    0.2,             0,                        mw,      -10 ,      -0.1*maxAnkle,        -20; % Min
                   0.25 ,    2.5,      maxAnkle,                        Mw,       10 ,       0.1*maxAnkle,         20}; % Max
               disp('Tau_retio = 12');

    case '6N_tagaLike_2Ank_torques_symm_feedback_eq' %removed the higher order adaptation, equalized number of parameters (10)

        Mamp = [maxAnkle,...
            maxHip];
        mamp = 0*Mamp;
        
        mw = 0*ones(1,4);
        Mw = maxW*ones(1,4);
         
        Keys = {'\tau_r','\tau_ratio', 'beta', 'amp_6n_semisymm',  '6neuron_taga_like_symm', 'k_hip_fb', 'IC_matsuoka';
                      1 ,       1,      1,             2,                         4,         1,            0 };
        Range = {  0.02 ,     0.1,    0.2,          mamp,                        mw,       -20; % Min
                   0.25 ,      50,    2.5,          Mamp,                        Mw,        20}; % Max
            
               disp('Tau_retio = 12');
       
    case '6N_tagaLike_2Ank_torques_symm_feedback_eq_adaptation' % need to update this

        Mamp = [maxAnkle,...
            maxHip];
        mamp = 0*Mamp;
        
        Madap = [1,1,100,100];
        madap = -Madap;
        
        mw = 0*ones(1,4);
        Mw = maxW*ones(1,4);
        
%         Keys = {'\tau_r','\tau_ratio', 'beta', 'amp_6n_semisymm',  '6neuron_taga_like_symm', 'k_hip_fb', 'k_adap' , 'dt_nom' , 'IC_matsuoka';
%                       1 ,       1,      1,             2,                         4,             1,          4     ,  1   0 };
%         Range = {  0.02 ,     0.1,    0.2,          mamp,                        mw,            -20,         madap,   0 ; % Min
%                    0.25 ,      50,    2.5,          Mamp,                        Mw,             20,         Madap,    1 }; % Max           
       Keys = {'\tau_r','\tau_ratio', 'beta', 'amp_6n_semisymm',  '6neuron_taga_like_symm', 'k_hip_fb', 'k_adap' , 'IC_matsuoka';
                      1 ,       1,      1,             2,                         4,             1,          4   ,  0 };
        Range = {  0.02 ,     0.1,    0.2,          mamp,                        mw,            -20,         madap; % Min
                   0.25 ,      50,    2.5,          Mamp,                        Mw,             20,         Madap}; % Max           

                   disp('Tau_retio = 12');

    case '6N_general_2Ank_torques'

        Mamp = [maxAnkle*ones(1,2*nAnkle1),...
            maxAnkle*ones(1,2*nAnkle2),...
            maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        mw = 0*ones(1,30);
        Mw = maxW*ones(1,30);
        
        Keys = {'\tau_r', 'beta', 'amp_6n_symm',                 'weights', 'ks_\tau',     'ks_c_6n_symm', 'IC_matsuoka';
                      1 ,      1,             1,                        30,        1 ,                 1 ,            0 };
        Range = {  0.02 ,    0.2,             0,                        mw,      -10 ,      -0.1*maxAnkle; % Min
                   0.25 ,    2.5,      maxAnkle,                        Mw,       10 ,       0.1*maxAnkle}; % Max
                      disp('Tau_retio = 12');
        
    case '6N_general_2Ank_torques_feedback'

        Mamp = [maxAnkle*ones(1,2*nAnkle1),...
            maxAnkle*ones(1,2*nAnkle2),...
            maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        mw = 0*ones(1,30);
        Mw = maxW*ones(1,30);
        
        % % Same structure as the previuos 6N CPG
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        Keys = {'\tau_r', 'beta', 'amp_6n_symm',                 'weights', 'ks_\tau',     'ks_c_6n_symm','k_hip_fb', 'IC_matsuoka';
                      1 ,      1,             1,                        30,        1 ,                 1 ,         1,         0 };
        Range = {  0.02 ,    0.2,             0,                        mw,      -10 ,      -0.1*maxAnkle,        -20; % Min
                   0.25 ,    2.5,      maxAnkle,                        Mw,       10 ,       0.1*maxAnkle,         20}; % Max
                       disp('Tau_retio = 12');
     
%     case '2_level_CPG'   % old version, includes parameters for
%     % adaptation, tonic input parameter (meaningless here)
%         Mamp = [maxAnkle*ones(1,2*nAnkle1),...
%             maxAnkle*ones(1,2*nAnkle2),...
%             maxHip*ones(1,2*nHip)];
%         mamp = 0*Mamp;
%         
%         % % Same structure as the 2N CPG
%         % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
%         % Include amp, offset, and duration of rectangular pulses for the PG
%         Mw = maxW;
%         mw = 0*Mw;
%         
%         nP = 4;
%         
%         MA = 20*ones(1,nP);
%         mA = -MA;
%         MD = ones(1,nP);
%         mD = 0*MD;
%         
%         N = 1;
%         Keys = {'\tau_r', 'beta','amp_2n_same_inputs',    '2neuron_symm_weights', 'ks_\tau',       'k_hip_fb',     'ks_c_2n_symm',  'PGamp',  'PGoffset', 'PGduration', 'IC_2_lvl';
%                       1 ,      1,                   1,                1,                 1 ,         1,             1,               nP,       nP,          nP,            0 };
%         Range = {  0.02 ,    0.2,                  mw,                1.084,           -10,        -20,            -20,              mA,       mD,          mD; % Min
%                    0.25 ,    20,                   Mw,                10,               10 ,        20,             20,              MA,       MD,          MD}; % Max      
%                
    case '2_level_CPG' % not up to date
        Mamp = [maxAnkle*ones(1,2*nAnkle1),...
            maxAnkle*ones(1,2*nAnkle2),...
            maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        % % Same structure as the 2N CPG
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        % Include amp, offset, and duration of rectangular pulses for the PG
        
        nP = 4;
        
        Mamp = 20*ones(1,nP);
        mamp = -Mamp;
        MD = ones(1,nP);
        mD = 0*MD;
        
        N = 1;
        Keys = {'\tau_r', 'beta',  '2neuron_symm_weights',       'k_hip_fb',    'PGamp',  'PGoffset', 'PGduration', 'IC_2_lvl';
                      1 ,      1,             1,                  1,              nP,       nP,          nP,            0 };
        Range = {  0.02 ,    0.2,         1.084,            -20,              mamp,       mD,          mD; % Min
                   0.25 ,    20,             10,                20,              Mamp,       MD,          MD}; % Max      
                         disp('Tau_retio = 12');
     
               
    case '2_level_CPG_eq'
        
        % % Same structure as the 2N CPG
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        % Include amp, offset, and duration of rectangular pulses for the PG
        Mw = maxW;
        mw = 0*Mw;
        
        nPA = 1;
        nPH = 1;
        nP = nPA+nPH;
        
        Mamp = 20*ones(1,nP);
        mamp = -Mamp;
        MD = ones(1,2*nP);
        mD = 0*MD;
        
        mPA = repmat([mamp(1),mD(1:2)],1,nPA);
        MPA = repmat([Mamp(1),MD(1:2)],1,nPA);
        mPH = repmat([mamp(1),mD(1:2)],1,nPH);
        MPH = repmat([Mamp(1),MD(1:2)],1,nPH);
        
        N = 1;
        Keys = {'\tau_r', 'beta',  '2neuron_symm_weights',    'k_hip_fb',   'PulseAnk', 'PulseHip', 'IC_2_lvl';
                      1 ,      1,              1,                  1,          3*nPA,     3*nPH,  0};
        Range = {  0.02 ,    0.2,          1.084,                -20,          mPA,       mPH; % Min
                   0.25 ,    20,              10,                 20,          MPA,       MPH}; % Max      
                            disp('Tau_retio = 12');
  
               
    case '2_level_CPG2' % not up to date
        Mamp = [maxAnkle*ones(1,2*nAnkle1),...
            maxAnkle*ones(1,2*nAnkle2),...
            maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        % % Same structure as the 2N CPG
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        % Include amp, offset, and duration of rectangular pulses for the PG
        Mw = maxW;
        mw = 0*Mw;
        
        nP = 4;
        
        Mamp = 20*ones(1,nP);
        mamp = -Mamp;
        MD = ones(1,nP);
        mD = 0*MD;
        
        N = 1;
        Keys = {'\tau_r', 'beta','amp_2n_same_inputs',    '2neuron_symm_weights', 'ks_\tau',       'k_hip_fb',     'ks_c_2n_symm',  'PGamp',  'PGoffset', 'PGduration', 'IC_2_lvl';
                      1 ,      1,                   1,                1,                 1 ,         1,             1,               nP,       nP,          nP,            0 };
        Range = {  0.02 ,    0.2,                  mw,                1.084,           -10,        -20,            -20,              mamp,       mD,          mD; % Min
                   0.25 ,    20,                   Mw,                10,               10 ,        20,             20,              Mamp,       MD,          MD}; % Max      
                            disp('Tau_retio = 12');
  
    case '2_level_CPG3' % not up to date
        
        Mamp = [maxAnkle*ones(1,2*nAnkle1),...
            maxAnkle*ones(1,2*nAnkle2),...
            maxHip*ones(1,2*nHip)];
        mamp = 0*Mamp;
        
        % % Same structure as the 2N CPG
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio) 
        % Include amp, offset, and duration of rectangular pulses for the
        % PG and frequency adaptation parameter k_om
        Mw = maxW;
        mw = 0*Mw;
        
        nP = 4;
        
        Mamp = 20*ones(1,nP);
        mamp = -Mamp;
        MD = ones(1,nP);
        mD = 0*MD;
        
        N = 1;
        Keys = {'\tau_r', 'beta','amp_2n_same_inputs',    '2neuron_symm_weights', 'ks_\tau',       'k_hip_fb',     'ks_c_2n_symm','k_om' ,  'PGamp',  'PGoffset', 'PGduration', 'IC_2_lvl';
                      1 ,      1,                   1,                1,                 1 ,         1,             1,              1         nP,       nP,          nP,            0 };
        Range = {  0.02 ,    0.2,                  mw,                1.084,           -10,        -20            -20               0         mamp,       mD,          mD; % Min
                   0.25 ,    20,                   Mw,                10,               10 ,        20             20               1         Mamp,       MD,          MD}; % Max      
           disp('Tau_retio = 12');

    case 'ConSpitz'

        % Include amp, offset, and duration of rectangular pulses for the PG
        Mw = maxW;
        mw = 0*Mw;
        
        nP = 4;
        
        Mamp = maxHip*ones(1,nP);
        mamp = -Mamp;
        MD = ones(1,nP);
        mD = 0*MD;
        
        N = 1;
        Keys = { 'omega' 'PGamp',  'PGoffset', 'PGduration', 'IC_1_lvl';
                      1      nP,       nP,          nP,            0 };
        Range = {   0.5      mamp,       mD,          mD; % Min
                      2      Mamp,       MD,          MD}; % Max
        
    case 'ConSpitz_eq'
 
        % Include amp, offset, and duration of rectangular pulses for the PG
        Mw = maxW;
        mw = 0*Mw;
   
        nPA = 1;
        nPH = 2;
        nP = nPA+nPH;
        
        Mamp = 20*ones(1,nP);
        mamp = -Mamp;
        MD = ones(1,2*nP);
        mD = 0*MD;
        
        mPA = repmat([mamp(1),mD(1:2)],1,nPA);
        MPA = repmat([Mamp(1),MD(1:2)],1,nPA);
        mPH = repmat([mamp(1),mD(1:2)],1,nPH);
        MPH = repmat([Mamp(1),MD(1:2)],1,nPH);
        
        N = 1;
        Keys = { 'omega',   'PulseAnk', 'PulseHip', 'IC_1_lvl';
                      1        3*nPA,     3*nPH,         0 };
        Range = {   0.5          mPA,       mPH; % Min
                      2          MPA,       MPH}; % Max
                  
    case 'ConSpitz_eq_adaptation'
 
        % Include amp, offset, and duration of rectangular pulses for the PG
        Mw = maxW;
        mw = 0*Mw;
   
        nPA = 1;
        nPH = 2;
        nP = nPA+nPH;
        
        Mamp = 20*ones(1,nP);
        mamp = -Mamp;
        MD = ones(1,2*nP);
        mD = 0*MD;
        
        Madap = [1,1,100,100];
        madap = -Madap;
        
        mPA = repmat([mamp(1),mD(1:2)],1,nPA);
        MPA = repmat([Mamp(1),MD(1:2)],1,nPA);
        mPH = repmat([mamp(1),mD(1:2)],1,nPH);
        MPH = repmat([Mamp(1),MD(1:2)],1,nPH);
        
        N = 1;
        Keys = { 'omega',   'PulseAnk', 'PulseHip','k_adap', 'IC_1_lvl';
                      1        3*nPA,     3*nPH, 4,         0 };
        Range = {   0.5          mPA,       mPH, madap; % Min
                      2          MPA,       MPH, Madap}; % Max
    case 'ConSpitz2'

        
        % % Same structure as the 2N CPG
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio)
        % Include amp, offset, and duration of rectangular pulses for the PG
        Mw = maxW;
        mw = 0*Mw;
        
        nP = 4;
        
        Mamp = 20*ones(1,nP);
        mamp = -Mamp;
        MD = ones(1,nP);
        mD = 0*MD;
        
        N = 1;
        Keys = { 'PGamp',  'PGoffset', 'PGduration', 'IC_1_lvl';
                      nP,       nP,          nP,            0 };
        Range = {     mamp,       mD,          mD; % Min
                      Mamp,       MD,          MD}; % Max
    case 'ConSpitz3'

        
        % % Same structure as the 2N CPG
        % Final genome with tau_r + beta (constant tau_u/tau_v ratio)
        % Include amp, offset, and duration of rectangular pulses for the PG
        Mw = maxW;
        mw = 0*Mw;
        
        nP = 4;
        
        Mamp = 20*ones(1,nP);
        mamp = -Mamp;
        MD = ones(1,nP);
        mD = 0*MD;
        
        N = 1;
        Keys = { 'k_om' 'PGamp',  'PGoffset', 'PGduration', 'IC_1_lvl';
                      1      nP,       nP,          nP,            0 };
        Range = {     0      mamp,       mD,          mD; % Min
                      1      Mamp,       MD,          MD}; % Max       
    otherwise
        error('invalid input');
end


for i=1:(length(Keys)-1)
    disp([Keys{1,i},' :']);
    disp(['    number of elements: ',num2str(Keys{2,i})]);
    disp(['    min: [',num2str(Range{1,i}),']']);
    disp(['    Max: [',num2str(Range{2,i}),']']);
end

save(genome_file, 'nAnkle1','nAnkle2', 'nHip', 'maxAnkle', 'maxHip', ...
    'Mamp', 'mamp', 'N', 'Mw', 'mw', ...
    'MutDelta0', 'MutDelta1', 'Keys', 'Range');

end

