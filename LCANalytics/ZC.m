function [eventvalue,isterminal,eventdir] = ZC(~,x,a,b,c,tau,T)

te=[x(1),x(1),x(2),x(2),x(1)-x(2)];

eventvalue=te;
isterminal=zeros(1,length(te));
eventdir=[1,-1,1,-1,-1];
end