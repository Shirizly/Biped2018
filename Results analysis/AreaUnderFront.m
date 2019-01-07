function A = AreaUnderFront(Fit,varargin)

%% receives Fit - matrix of first front fitness values over 2 fitness functions
%% scale - 2d vector of scaling proportions for the 2 fitnesses.
%% returns normalized area under the front
% Fit =    [0.1 4;...
%     0.4 3.6;...
%     0.5 3.2;...
%     0.6 2.6;...
%     0.75 2.1;...
%     0.8 1.6;...
%     0.85 1.4;...
%     0.9 1.3;...
%     0.98 1.1];

% Scale = [1 12];

% if nargin<2
%     Scale = max(Fit);
% else
%     Scale = varargin{1};
% end
% normalv = repmat(Scale,size(Fit,1),1);
NFit = Fit;%./normalv;
NFit = sortrows(NFit);



if any(any(isnan(NFit)))||size(NFit,1)<2
    A = 0;
else
    mf1 = min(NFit(:,1));
    mf2 = min(NFit(:,2));
    if mf1 >0
        Mf2 = NFit(find(NFit(:,1)==mf1,1,'first'),2);
        NFit = [0,Mf2 ; NFit];
    end
    
    Mf1 = max(NFit(:,1));
    Mf2 = max(NFit(:,2));
    
    
    X = NFit(:,1);
    Y = NFit(:,2);
    
    try
    f1 = fit(X,Y,'linearinterp');
    
    %         plot(X,Y)
    
        A = integrate(f1,Mf1,0);
        A = A/(Mf1*Mf2);
    catch err
        blip = 1;
    end
end



end


% Fit =    [1.6904         0         0;...
%     1.5640    0.0277    0.0170;...
%     1.5111    0.0472    0.1169;...
%     1.5093    0.0365    1.1231;...
%     1.4565    0.0410    0.3060;...
%     1.4498    0.0528    1.0163;...
%     1.4110    0.0551    0.1934;...
%     1.3888    0.0556    0.2523;...
%     1.3769    0.0575    0.1319];
% Fit = Fit(:,2:3);
% scale = [1 10];