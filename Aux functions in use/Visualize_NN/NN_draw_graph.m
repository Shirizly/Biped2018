function NN_draw_graph( net,inputsNames)
% this function takes a NN and plot a color maps of its weights

% IMPORTANT: this function is for feedfarward NN only!

% inputs:
% *) 'net' - the NN structure
% *) 'inputNames' - the names of the features in the input array X
%                   reminder- the features are the rows, the samples are
%                   the colums.
% *) 'graphType' - options:
%                   #) 'grey scale abs' - show the weightMap in grey scale
%                               the larger abs(W_ij) the darker the square.
%                   #) 'red and blue' - W_ij>0 and increasing the more red
%                               it is. if W_ij<0 and decreasing more blue.
% 
% 
% 
% 
% outputs:
% none

% Some more information:
% *) 'net.IW' - 'numLayers'-by-'numInputs' cell array of input weight values
%         note: 'numInputs' - is note the nuber of features in Xin, 
%                     but rather the number of vectors Xin (usually = 1)
%         each cell in 'IW' represents:
%         {[vec1];   = from the inputs to the first layer
%          [vec2];   = from the input to the 2nd layer...
%          ...}
% *) 'net.LW' - 'numLayers'-by-'numLayers' cell array of layer weight values
%                 LW_ij represents connection from layer 'j' to layer 'i'
%         each cell in 'LW' represents:
%         {[1st layer to itself], [2nd layer to the 1st];
%          [1st layer to the 2nd], [2nd layer to itself]};...
% *) 'net.b' - 'numLayers'-by-1 cell array of bias values
%             represents the bias of each layer

% Net Parameters:
inputsNum = size(net.IW{1,1},2); % number of inputs
layerNum = numel(net.b); % number of layers                  
nodesEdges_Xcoor = cell(1,layerNum+1);
nodesEdges_Ycoor = cell(1,layerNum+1);

% Plot Parameters:
Xmin = 0; % x min for the plotting area;
Xmax = 1; % x max for the plotting area;
Ymin = 0;
Ymax = 1;

% get the NN arcitecture:
architecture = zeros(1,layerNum);
for l=1:layerNum
    architecture(1,l) = size(net.b{l,1},1);
end

% declare figure
figure; hold on;

% Draw the inputs:
inputNodesR = 0.15 * (Ymax - Ymin)/(inputsNum+2); 
inputNodes_Xcoor = Xmin+inputNodesR;
inputNodes_Ycoor = linspace(Ymin,Ymax,inputsNum+2);
for in=1:inputsNum
    % drawNode(xCoor,yCoor,Radius)
    drawNode(inputNodes_Xcoor,...
        inputNodes_Ycoor(1,in+1),...
        inputNodesR)
end
% save the coor. of the node edge:
nodesEdges_Ycoor{1,1} = inputNodes_Ycoor(1,2:end-1);
nodesEdges_Xcoor{1,1} = inputNodes_Xcoor;

% Draw the neurons nodes:
NeuronNodesR = 0.3 * (Ymax - Ymin)/(max(architecture)+2);
NeuronNodes_Xcoor = linspace(Xmin,Xmax,(length(architecture)+2));
for l=1:layerNum
    NeuronNodes_Ycoor = linspace(Ymin,Ymax,architecture(1,l)+2);
    for n=1:architecture(1,l)
        % drawNode(xCoor,yCoor,Radius)
        drawNode(NeuronNodes_Xcoor(1,l+1),...
            NeuronNodes_Ycoor(1,n+1),...
            NeuronNodesR)
    end
    nodesEdges_Ycoor{1,l+1} = NeuronNodes_Ycoor(1,2:end-1);
    nodesEdges_Xcoor{1,l+1} = NeuronNodes_Xcoor(1,l+1);
end


% Draw the weights:
for l=1:layerNum
    if l==1
        weights = net.IW{1,1};
        X = [nodesEdges_Xcoor{1,l}+inputNodesR,... % line origin point
            nodesEdges_Xcoor{1,l+1}-NeuronNodesR];
    else
        weights = net.LW{l,l-1};
        X = [nodesEdges_Xcoor{1,l}+NeuronNodesR,...
            nodesEdges_Xcoor{1,l+1}-NeuronNodesR];
    end
    
    % plot weights as lines:
    for n1 = 1:size(nodesEdges_Ycoor{1,l},2)
        % run on the outputs of previuos layer (note: if 'l=1' then this is the NN inputs)
        for n2 = 1:size(nodesEdges_Ycoor{1,l+1},2)
            % run on the nerons in the current layer
            
            % set 'LineWidth' as the abs(weight):
            width = abs(weights(n2,n1));
            % set color as the value of W_ij (blue- inhibitory,
            % red-excetatory, wight- zero weigth) :
            if weights(n2,n1)>0
                COLOR = 'r';
            else
                COLOR = 'b';
            end
            Y = [nodesEdges_Ycoor{1,l}(1,n1),...
                nodesEdges_Ycoor{1,l+1}(1,n2)];
            drawWeigth(X,Y,width,COLOR)
        end
    end
end

% set layer names as Xaxis Ticks:
set(gca,'XTick',cell2mat(nodesEdges_Xcoor));
for l=1:layerNum+1
    if l==1
        LABELS{1,l} = 'inputs';
    elseif l==layerNum+1
        LABELS{1,l} = 'outputs';
    else
        LABELS{1,l} = sprintf('hidden_%d',l-1);
    end
end
set(gca,'XTickLabel',LABELS);

% set inputs names as Yaxis Ticks:
set(gca,'YTick',nodesEdges_Ycoor{1,1});
set(gca,'YTickLabel',inputsNames);

title({'Network architecture:',...
    'blue = negative weight   ;   red = positive weight'});

    function drawNode(x,y,R)
        % 'CircRes' - resolution of the cericle drawing
        CircRes = 10;
        
        coordX=zeros(1,CircRes);
        coordY=zeros(1,CircRes);
        coordZ=zeros(1,CircRes);
        
        for r=1:CircRes
            coordX(1,r) = x + R*cos(r/CircRes*2*pi);
            coordY(1,r) = y + R*sin(r/CircRes*2*pi);
            coordZ(1,r) = 0;
        end
        h = patch(coordX,coordY,coordZ,'k');
        set(h,'EdgeColor',[0,0,0]);
        set(h,'LineWidth',1);
    end

    % Draws a vector from x0 to x1
    function drawWeigth(x0,x1,Lwidth,Lcolor)
        % 'x0' - origin point
        % 'x1' - end point
        plot(x0,x1,'LineWidth',Lwidth,'Color',Lcolor);
    end
end

