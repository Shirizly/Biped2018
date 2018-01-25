function [Fronts] = Pareto(GA, Data, Inv, NFrontsReq) %#ok<INUSL>
%PARETO Finds Pareto fronts for multi-objective genetic algorithms
%   The algorithm builds each front by "dropping" elements that are
%   dominated until only dominant or weak dominant elements remain.
%   These elements are removed from the set and the process is repeated to
%   find a second pareto front and so forth until all elements are
%   accounted for.
%   When Inv = 1, domination is defined by >= instead of <=.
%   NOTE: Data should be provided a unique ID as the last column.

% Version 1.1 - 28/08/2014
if nargin<4
    NFrontsReq = [];
end

if nargin<3
    Inv = 0;
end

if nargin<2
    nI = 5; nJ = 1000;
    % Run sample code
    % Build sample data
    Data = zeros(nI*nJ,3);
    for i = 1:nI
%         r = 10*i;
        for j = 1:nJ
%             r = 10*i+rand();
            r = randsample(2:2:10,1);
            phi = -pi/2+pi*rand();
            th = 2*pi*rand();
            Data(nJ*(i-1)+j,:) = [...
                r*cos(phi)*cos(th), ...
                r*cos(phi)*sin(th), ...
                r*sin(phi)];
        end
    end
    Data = abs(Data);
    
    % Plot sample data before running algorithm
%     figure()
%     plot3(Data(:,1),Data(:,2),Data(:,3),'+');
    
    % Save data to origx
    origx = Data;
    
    % Give unique ID to each sample
    Data = [Data (1:size(Data,1))'];
end

% Round off to 3 decimal places
% x = round(x*1000)/1000;

% Sort by first objective
if Inv
    SortedData = sortrows(Data,1);
else
    SortedData = sortrows(Data,-1);
end

if size(Data,2) == 1
    Fronts = num2cell(SortedData(:,end));
    return;
end

% Start separating layers
NFronts = 1;
Fronts = {};
out = [];

while 1
    % Drop samples that are "dominated"
    % i.e. that x1>x2 for each objective
    i = 1;
    Nd = size(SortedData,1);
    if Nd == 0
        break;
    end
    
    while i < Nd
        j = i+1;
        
        % Find dominated values efficiently :)
        if Inv
            IDs = j-1 + find(all(SortedData(j:Nd,2:end-1) >= ...
                            repmat(SortedData(i,2:end-1),Nd-j+1,1),2)==1);
        else
            IDs = j-1 + find(all(SortedData(j:Nd,2:end-1) <= ...
                            repmat(SortedData(i,2:end-1),Nd-j+1,1),2)==1);
        end
        SortedData(IDs,:) = [];
        Nd = size(SortedData,1);
        
        i = i+1;
    end
    
    % Save "undominated" samples to current front
    Fronts{NFronts} = SortedData(:,end); %#ok<AGROW>
    out = [out; Fronts{NFronts}]; %#ok<AGROW>
    NFronts = NFronts+1;
    
    if ~isempty(NFrontsReq) && NFronts>NFrontsReq
        break;
    end
    
    % Restore samples to original data minus samples already on fronts
    SortedData = Data; SortedData(ismember(SortedData(:,end),out),:) = [];
    if Inv
        SortedData = sortrows(SortedData,1);
    else
        SortedData = sortrows(SortedData,-1);
    end
end

if nargin<2
    % Run sample code
    close all
    N=length(Fronts);

    figure()
    hold on
    for i=1:N
        if length(Fronts{i})<3
            continue;
        end
        Data = origx(Fronts{i},1);
        y = origx(Fronts{i},2);
        z = origx(Fronts{i},3);
        tri = delaunay(Data,y);
        trisurf(tri, Data, y, z,'FaceColor',[1-i/N; i/2/N; i/N]);
    end
end

end

%% %%%%%%%%%%%%%% Slow method %%%%%%%%%%%%%% %%
% % Get number of samples and objectives
% data = size(x,1);
% % Give unique ID to each sample
% x=[x (1:data)'];
% tempx = x;
% NFronts=1;
% Fronts={};
% out = [];
% while 1
%     i = 1;
%     while i <= data
%         j = 1;
%         while j <= data
%             if i == j
%                 j = j+1;
%                 continue;
%             end
% 
%             if tempx(j,1:end-1) <= tempx(i,1:end-1)
%                 tempx(j,:) = [];
%                 data = data - 1;
%                 if i>j
%                     i = i-1;
%                 end
%             else
%                 j = j+1;
%             end
%         end
%         i = i+1;
%     end
% 
%     Fronts{NFronts} = tempx(:,end); %#ok<AGROW>
%     out = [out; Fronts{NFronts}]; %#ok<AGROW>
%     NFronts = NFronts+1;
% 
%     tempx = x; tempx(out,:) = [];
%     data = size(tempx,1);
%     if data == 0
%         break;
%     end
% end

