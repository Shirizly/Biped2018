function NN_weights_matrix_plot( net,whichOne,inputNames,graphType )
% this function takes a NN and plot a color maps of its weights

% IMPORTANT: this can be replace with "heatmap" function!!

% inputs:
% *) 'net' - the NN structure
% *) 'whichOne' - either:
%                 {'IW',i} - weights from Input to layer 'i'
%                 {'LW',i,j} - weights from Layer 'j' to Layer 'i'
%                 {'IW',i,'b'} or {'LW',i,j,'b'} -  also show bias
%       Example: enter cell array like:
%               #) {'IW',1} - show weights from Input to the 1st Layer
%               #) {'LW',3,2} - show weights from Layer 2 to the Layer 3
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

% % get the NN weights
switch whichOne{1,1}
    case 'IW'   % show weights from input to layer 'i'
        weights = net.IW{whichOne{1,2},1};
        
    case 'LW'   % {'LW',i,j} - weights from Layer 'j' to Layer 'i'
        weights = net.LW{whichOne{1,2},whichOne{1,3}};
        
    otherwise
        error('unknown input...');
end

% % Make weights names and add sub labels and figure title
switch whichOne{1,1}
    case 'IW' % names are the names of the inputs
        parametersCells = inputNames;
        Text_under_Xlabel = 'Names of the net inputs';
        Text_nextTo_Ylabel = sprintf('To neuron_i at layer %d',whichOne{1,2});
        Title = ...
            sprintf('NN weights from the input to layer %d',whichOne{1,2});
        
    case {'LW','b'} % names contain the Layer index
        Text_under_Xlabel = sprintf('From neuron_j at layer %d',whichOne{1,3});
        Text_nextTo_Ylabel = sprintf('To neuron_i at layer %d',whichOne{1,2});
        Title = ...
            sprintf('NN weights from layer %d to layer %d',...
            whichOne{1,3},whichOne{1,2});
        
        parametersCells = cell(1,size(weights,2));
        for n=1:size(weights,2)
            parametersCells{1,n} = ...
                sprintf('y^{(%d)}_{%d}',whichOne{1,3},n);
            % 'y ^l _j'-output of neuron 'j' at Layer 'l'
        end
end

% % if we want to also show the layer bias:
if strcmp(whichOne{1,end},'b')
    bias = net.b{whichOne{1,2},1};
    weights = [weights,bias]; % add the bias vector to the weights
    parametersCells{1,end+1} = sprintf('layer bias');
end

% % make neurons names: (in the recieving layer)
for n=1:size(weights,1)
    HiddenNumCells{1,n} = ...
        sprintf('N^{(%d)}_{%d}',whichOne{1,2},n);
    % 'N ^l _j'- nueron num 'j' at Layer 'l'
end

[x,y] = meshgrid(1:size(weights,2),1:size(weights,1));   %# Create x and y coordinates for the strings
textStrings = num2str(weights(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding

figure;
switch graphType
    case 'grey scale abs'
        imagesc(abs(weights));
        colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                                 %#   black and lower values are white)
        hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
        midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
        textColors = repmat(abs(weights(:)) > midValue,1,3);%# Choose white or black for the
                                                     %#   text color of the strings so
                                                     %#   they can be easily seen over
                                                     %#   the background color
        set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors

    case 'red and blue'
        % creating diverging color map:
        mapResolution = 100; % how many colors
        [MAP] = diverging_map(0:1/mapResolution:1,...
            [40,53,147]./255,[145,28,8]./255);
        colormap(MAP);
        
        M = abs(max(max(weights)));
        m = abs(min(min(weights)));
        if M-m>=0      
           clim=[-m m];
           else
           clim=[-M M];
        end
        
        imagesc(weights,[clim]); colorbar;

        hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
        midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
        
        % leave texts on light colors black. but set texts on light colors
        % white:
        % find cells with light backround colors by checking if the value
        % is in a range next to the mean:
        Cdelta = 1; % set the range limits
        lightColors_ind = (weights(:) > midValue+Cdelta) |...
            (weights(:) < midValue-Cdelta);
        
        textColors =...
            repmat(lightColors_ind,1,3);%# Choose white or black for the
                                                     %#   text color of the strings so
                                                     %#   they can be easily seen over
                                                     %#   the background color
        set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
end

title(Title);
set(gca,'XTick',1:size(weights,2),...     %# Change the axes tick marks
        'XTickLabel',parametersCells,...  %#   and tick labels
        'YTick',1:size(weights,1),...
        'YTickLabel',HiddenNumCells,...
        'TickLength',[0 0]);
% set additional text under the X-label:    
text(0.5, -0.1, Text_under_Xlabel,...
    'Units', 'normalized','HorizontalAlignment','center');
% set additional text under the X-label:    
htext=text(-0.1, 0.5, Text_nextTo_Ylabel,...
    'Units', 'normalized','HorizontalAlignment','center');
set(htext,'Rotation',90);


end

