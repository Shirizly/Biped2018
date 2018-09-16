function [ppv,dppv] = Sto_Ter_Array_Gen(FileName,i,runs,TerVarS,TerVarE,nseg,xend)
%% Create Sto-Terrains for GA runs:

% runs = 10
% 
% 
% slo_num = 10
% Slo_lim = 0.3
% xend = 10


pps = cell(runs);
dpps = cell(runs);
for j = 1:runs
%   [ppv{j} dppv{j}] = sto_ter_gen(slo_num,Slo_lim,xend);
%     [ppv{j} dppv{j}] = Sto_Ter_Gen2(TerVar,xend);
    [ppv{j} dppv{j}] = Sto_Ter_Gen4(TerVarS,TerVarE,nseg,xend);
end

save([FileName num2str(i)],'ppv','dppv','TerVarS','TerVarE','nseg','xend');
end