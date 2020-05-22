function Population = EnvironmentalSelection(Population,N)
    %Update Rank & select N fittness member
    PopObj = Population.objs;
    NF  = NDSort(PopObj,inf);
    CD = -CrowdingDistance(PopObj,NF);
    [~,Rank] = sortrows(cat(2,NF',CD'));
    Population = Population(Rank(1:N));
    for i = 1:N
        Population(i).rank=i;
    end

end