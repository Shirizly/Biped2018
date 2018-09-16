function [pmin,pmax] = Find_limits(data)

% Start separating layers

bound = boundary(data(:,1),data(:,2),0); %

datawork = (data(bound,:));

ext = [max(datawork(:,1)) min(datawork(:,1))];

p(1) = find(datawork(1:end-1,1) == ext(1),1,'first');
p(2) = find(datawork(1:end-1,1) == ext(2),1,'first');
    
p = sort(p);
b1 = datawork(p(1):p(2),:);
b2 = datawork([p(2):end-1,1:p(1)],:);

flag = 1;
while flag
    flag = 0;
    if b1(1,1)==b1(2,1)
        b1(1,:)=[];
        flag = 1;
    end
    if b1(end,1)==b1(end-1,1)
        b1(end,:)=[];
        flag = 1;
    end
    if b2(1,1)==b2(2,1)
        b2(1,:)=[];
        flag = 1;
    end
    if b2(end,1)==b2(end-1,1)
        b2(end,:)=[];
        flag = 1;
    end
end

if mean(b1(:,2))>mean(b2(:,2))
    pmax = sortrows(b1);
    pmin = sortrows(b2);
else
    pmax = sortrows(b2);
    pmin = sortrows(b1);
end
% figure
% hold on
% plot(data(:,1),data(:,2),'+');
% plot(datawork(:,1),datawork(:,2));
% hold off





end

% Nd = size(datawork,1);
% if Nd == 0
%     return;
% end
% i = 1;
% 
% while i < Nd
%     j = i+1;
%     IDs = j-1 + find(all(datawork(j:Nd,2:end) >= ...
%         repmat(datawork(i,2:end),Nd-j+1,1),2)==1);
%     datawork(IDs,:) = [];
%     Nd = size(datawork,1);
%     
%     i = i+1;
% end
% 
% 
% pmin = datawork;
% 
% datawork = sortrows(data(bound,:),-1);
% 
% figure
% hold on
% plot(data(:,1),data(:,2),'+');
% plot(datawork(:,1),datawork(:,2));
% hold off
% 
% 
% Nd = size(datawork,1);
% if Nd == 0
%     return;
% end
% i = 1;
% 
% while i < Nd
%     j = i+1;
%     
%     % Find dominated values efficiently :)
%     IDs = j-1 + find(all(datawork(j:Nd,2:end) <= ...
%         repmat(datawork(i,2:end),Nd-j+1,1),2)==1);
%     
%     datawork(IDs,:) = [];
%     Nd = size(datawork,1);
%     
%     i = i+1;
% end
% 
% 
% pmax = datawork;
