function Population = subNSGAII(Population,Tasks,rmp)
    Global  = GLOBAL.GetObj();
    while Global.NotTermination(Population)
        FrontNo  = NDSort(Population.objs,inf);
        CrowdDis = CrowdingDistance(Population.objs,FrontNo);
        ParentPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = OffspringCreation(Population(ParentPool),Tasks,rmp);
        Population = EnvironmentalSelection([Population,Offspring],Global.N);
        if mod(Global.gen,10) == 0 
            Population = UpdateTask(Population,Tasks); 
        end 
    end  
end

