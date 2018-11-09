function ParGra(nI,nD,nJ) %#ok<INUSL>
%PARETO Finds Pareto fronts for multi-objective genetic algorithms
%   The algorithm builds each front by "dropping" elements that are
%   dominated until only dominant or weak dominant elements remain.
%   These elements are removed from the set and the process is repeated to
%   find a second pareto front and so forth until all elements are
%   accounted for.
%   When Inv = 1, domination is defined by >= instead of <=.
%   NOTE: Data should be provided a unique ID as the last column.


    NFrontsReq = [];

    Inv = 0;

%     nI = 5; nJ = 1000;
    % Run sample code
    % Build sample data
    Data = zeros(nI*nJ,nD);
    for i = 1:nI
%         r = 10*i;
        for j = 1:nJ
%             r = 10*i+rand();
            r = randsample(2:2:10,1)*(1+(rand-0.5)/20);
            phi = -pi/2+pi*rand();
            th = 2*pi*rand();
            if nD == 3
                Data(nJ*(i-1)+j,:) = [...
                    r*cos(phi)*cos(th), ...
                    r*cos(phi)*sin(th), ...
                    r*sin(phi)];
            else
                Data(nJ*(i-1)+j,:) = [...
                    r*cos(phi), ...
                    r*sin(phi)];
            end
        end
    end
    Data = abs(Data)/max(max(abs(Data)));

    
    
    % Plot sample data before running algorithm
%     figure()
%     plot3(Data(:,1),Data(:,2),Data(:,3),'+');
    
    % Save data to origx
    origx = Data;
    
    % Give unique ID to each sample
    Data = [Data (1:size(Data,1))'];



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


    % Run sample code
    close all
    N=length(Fronts)

    figure()
    hold on
    if nD == 3 
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
    else
        title('Pareto fronts of population','fontSize',20,'Interpreter','Latex')
        xlabel('1st Fitness Function','fontSize',16,'Interpreter','Latex')
        ylabel('2nd Fitness Function','fontSize',16,'Interpreter','Latex')
        for i=1:N
            if length(Fronts{i})<2
                continue;
            end
            Data = origx(Fronts{i},1);
            y = origx(Fronts{i},2);
            plot(Data,y,'-o','lineWidth',2,'markerSize',6,'Color',[i/N;i/2/N;1-i/N]);
        end
    end
    grid on
    hold off
end


% i/2/N


