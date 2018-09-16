function V = VolUnderFront(Fit,varargin)

%% receives Fit - matrix of first front fitness values over 2 fitness functions
%% scale - 2d vector of scaling proportions for the 2 fitnesses.
%% returns normalized area under the front
Fit(:,1) = max(Fit(:,1)-ones(size(Fit(:,1))),0);
if nargin<2
    Scale = max(Fit);
else
    Scale = varargin{1};
end
normalv = repmat(Scale,size(Fit,1),1);
NFit = Fit./normalv;
NFit = sortrows(NFit);
ord = [1 2 3;2 3 1;3 1 2];
Vol = zeros(3,1);


if any(any(isnan(NFit)))||size(NFit,1)<3
    V = 0;
else
   
    for i = 1:3
        
        X = NFit(:,ord(i,1));
        Y = NFit(:,ord(i,2));
        Z = NFit(:,ord(i,3));

        f1 = fit([X,Y],Z,'linearinterp');
        min1 = min(X);
        max1 = max(X);
        
        [pmin,pmax] = Find_limits([X,Y]);
        
        min2 = fit(pmin(:,1),pmin(:,2),'linearinterp');
        max2 = fit(pmax(:,1),pmax(:,2),'linearinterp');
        
        Vol(i) = quad2d(f1,min1,max1,min2,max2,'absTol',1e-3,'relTol',1e-1);
%         if Vol(i)>0.2
%             figure
%             plot(f1,[X,Y],Z);
%         end
    end
    V = mean(Vol);
end

% plot(f1,[X,Y],Z);

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
% scale = [1.1 1 10];