f = 1;
            IDsOut2 = IDabove{groupsd+1};
            while 1
                try
                FL = length(Fronts{f});
                catch me
                    save('error.mat',me);
                end
                if length(IDsOut2)-FL<Nkeep
                    break
                end
                IDsOut2 = setdiff(IDsOut2, Fronts{f});
                f = f+1;
            end
            