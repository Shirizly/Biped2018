function [ Result, Seq ] = CheckGenome( Ge, Seq )
% Checks a gene sequence to verify that it satisfies
% all the min/max conditions
    Result = {1,''};
    
    % Check that it's the correct length
    if length(Seq)~=Ge.Length
        Result = {0;
            ['Decode failed. Genes provided: ',num2str(length(Seq)),...
            '. Genes required: ',num2str(Ge.Length)]};
        return;
    end
    
    % Check that it's between the max and min range
    if any(Seq<Ge.Range(1,:) | Seq>Ge.Range(2,:))
        Result = {0;
        	'Sequence outside allowed genome range'};
        return;
    end
    beta = Seq(2);
    beta_pos = 2;

    SeqPos = 1;
    for k = 1:size(Ge.Keys,2)
        % Go over each key
        switch Ge.Keys{1,k}
            case 'Pulses'
                NPulses = Ge.Keys{2,k}(1);
                SeqPos0 = SeqPos;
                x = 0:0.001:1;
                Torque = zeros(1,length(x));
                for p = 1:NPulses
                    % Check that the pulse ends before the phase resets
                    End = Seq(SeqPos0+1)+Seq(SeqPos0+2);
                    if End>1
                        % Shorten both the offset and the duration
                        % so that End will be slightly smaller than 1
                        Seq(SeqPos0+1) = Seq(SeqPos0+1)/(End*1.001);
                        Seq(SeqPos0+2) = Seq(SeqPos0+2)/(End*1.001);
                        Result{2} = 'Pulse shortened';
                    end
                    
                    Torque = Torque + ...
                        Seq(SeqPos0)*(x>Seq(SeqPos0+1) & x<End);
                    SeqPos0 = SeqPos0 + Ge.KeyLength.Pulses;
                end
                % Check that the compound torque signal doesn't exeed
                % the min/max allowed
                factor = [];
                if any(Torque<Ge.Range(1,SeqPos)) % Min
                    factor = min(Torque)/Ge.Range(1,SeqPos);
                end
                if any(Torque>Ge.Range(2,SeqPos)) % Max
                    factor = max(Torque)/Ge.Range(2,SeqPos);
                end
                if ~isempty(factor)
                    % Fix the genome
                    SeqPos0 = SeqPos;
                    for p = 1:NPulses
                        Seq(SeqPos0) = Seq(SeqPos0)/factor;
                        SeqPos0 = SeqPos0 + Ge.KeyLength.Pulses;
                    end
                    Result{2} = ['Compound torque reduced by factor of ',num2str(factor)];
                end
                
            case {'PGamp','PGamp2','PGamp3'}
                NPulses = Ge.Keys{2,k};
                dp = Ge.Keys{2,k};
                for p = 1:NPulses
                    % Check that the pulse ends before the phase resets
                    End = Seq(SeqPos+dp+p-1)+Seq(SeqPos+p-1+dp*2);
                    if End>1
                        % Shorten the duration
                        % so that End will be slightly smaller than 1
                        Seq(SeqPos+p-1+dp*2) = 0.999 - Seq(SeqPos+dp+p-1);
                        Result{2} = 'Pulse shortened';
                    end                   
                end
                SeqPos = SeqPos + Ge.Keys{2,k}*2; % because decoding amp decode also offset and duration
                
            case {'PulseAnk','PulseHip'}
                seqlen = Ge.Keys{2,k};
                for p = 1:3:seqlen
                    % Check that the pulse ends before the phase resets
                    End = Seq(SeqPos+p)+Seq(SeqPos+p+1);
                    if End>1
                        % Shorten the duration
                        % so that End will be slightly smaller than 1
                        Seq(SeqPos+p+1) = 0.999 - Seq(SeqPos+p);
                        Result{2} = 'Pulse shortened';
                    end                   
                end
                if (seqlen/3)>1 % for 2 pulses case - ordering the pulses based on offset
                    if Seq(SeqPos+1)>Seq(SeqPos+4) % first pulse comes after scond pulse
                        tempseq = Seq(SeqPos:SeqPos+2); % reordering the sequence
                        Seq(SeqPos:SeqPos+2) = Seq(SeqPos+3:SeqPos+5);
                        Seq(SeqPos+3:SeqPos+5) = tempseq;
                    end
                end
                
            case 'ExtPulses'
                NPulses = Ge.Keys{2,k}(1);
                SeqPos0 = SeqPos;
                
                Total = 0;
                for p = 1:NPulses
                    Total = Total + Seq(SeqPos0);
                    SeqPos0 = SeqPos0 + Ge.KeyLength.ExtPulses;
                end
                
                factor = [];
                if Total<Ge.Range(1,SeqPos) % Min
                    factor = Total/Ge.Range(1,SeqPos);
                end
                if Total>Ge.Range(2,SeqPos) % Max
                    factor = Total/Ge.Range(2,SeqPos);
                end
                if ~isempty(factor)
                    % Fix the genome
                    SeqPos0 = SeqPos;
                    for p = 1:NPulses
                        Seq(SeqPos0) = Seq(SeqPos0)/factor;
                        SeqPos0 = SeqPos0 + Ge.KeyLength.ExtPulses;
                    end
                    Result{2} = ['Compound event triggered torque reduced',...
                        'by factor of ',num2str(factor)];
                end
                
            case '2neuron_symm_weights' % appears in the genetic sequence after beta, s relies on beta being read already
                a = Seq(SeqPos);
                if(a>beta+1)
                    Seq(beta_pos) = a - 1 + (0.01 + rand(1)*(Ge.Range(2,k) - a + 1));
                    Result{2} = 'beta parameter corrected to fulfill matsuoka conditions';
                end                
                T = tau*12;
                wn = 1/T*sqrt(((tau+T)*beta-tau*a)/(tau*a));
                wl = 2*pi;
                
                if wn > wl
                    Seq(tau_pos) = (wn/wl)*tau; %because the frequency depend on tau like wn ~ 1/tau
                end
                
            case 'beta' % keeping value and pos of beta for matsuoka condition check in case '2neuron_symm_weights'
                beta = Seq(SeqPos);
                beta_pos = SeqPos;     
            
            case '\tau_r'
                tau = Seq(SeqPos);
                tau_pos = SeqPos;
        end
        
        
        % Move the sequence reading position
        SeqPos = Ge.AdvSeq(SeqPos,k);
    end
end

