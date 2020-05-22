function Population = subMOEAD(Population,Tasks,rmp)
     %% Generate the weight vectors
    Global  = GLOBAL.GetObj();
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    T = ceil(Global.N/10);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
   
    Z = min(Population.objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        % For each solution
        for i = 1 : Global.N      
            % Choose the parents
            P = B(i,randperm(size(B,2)));

            % Generate an offspring
            Offsprings  = OffspringCreation(Population(P(1:2)),Tasks,rmp);
            Offspring = Offsprings(1);

            % Update the ideal point
            Z = min(Z,Offspring.obj);

            % Update the neighbours (PBI approach)
            normW   = sqrt(sum(W(P,:).^2,2));
            normP   = sqrt(sum((Population(P).objs-repmat(Z,T,1)).^2,2));
            normO   = sqrt(sum((Offspring.obj-Z).^2,2));
            CosineP = sum((Population(P).objs-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
            CosineO = sum(repmat(Offspring.obj-Z,T,1).*W(P,:),2)./normW./normO;
            g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
            g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
                
            Population(P(g_old>=g_new)) = Offspring;
        end
        Population = EnvironmentalSelection(Population,length(Population));% only update the rank
        if mod(Global.gen,10) == 0 
            Population = UpdateTask(Population,Tasks); 
        end
    end

end