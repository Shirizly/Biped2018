function CM = CellMap2Pi(x_range,y_range,z_range,map_fun)
nX = length(x_range);
nY = length(y_range);
nZ = length(z_range);
CM = -1*ones(nX-1,nY-1,nZ-1);
%%
ind0 = GetCell([0;0;0],x_range,y_range,z_range);
CM(ind0(1),ind0(2),ind0(3)) = 0;
% -1 - not assigned 
% >=0 - converge
% -2 - diverge

for i = 1:length(CM(:))
    if CM(i) ~= -1
        continue
    else
        ind = [i];
        stop = 0;
        [xx,yy,zz] = ind2sub(size(CM),ind(end));
        ic = PointFromCell([xx,yy,zz],x_range,y_range,z_range)';
        while ~stop
            ic = map_fun(ic);
            if sum(isinf(ic))==length(ic)
                    N=length(ind);
                    CM(ind) = linspace(N,1,N)*1i;
                    stop = 1;
            else
                cur_cell = GetCell(ic,x_range,y_range,z_range)
                    if ~isnan(cur_cell)
                        cur_ind = sub2ind(size(CM),cur_cell(1),cur_cell(2),cur_cell(3));
                        ind = [ind,cur_ind]; 
%                         display(cur_ind);
%                         display(CM(cur_ind));
%                         display(ic);
                        if CM(cur_ind) >= 0
                            N=length(ind);
                            CM(ind) = linspace(CM(cur_ind)+(N-1),CM(cur_ind),N);
                            stop = 1;
                        end
                        if imag(CM(cur_ind)) >0
                            N=length(ind);
                            CM(ind) = linspace(CM(cur_ind)+(N-1)*1i,CM(cur_ind),N);
                            stop = 1;
                        end
                    end
            end
        end
        disp([num2str(sum(CM(:)~=-1)),' out of ',num2str((nX-1)*(nY-1)*(nZ-1))]);
    end
end
       
            
            
                
    