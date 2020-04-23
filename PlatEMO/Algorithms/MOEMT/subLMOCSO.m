function Population = subLMOCSO(Population,Tasks,rmp)
    Global  = GLOBAL.GetObj();
    [V,Global.N] = UniformPoint(Global.N,Global.M);
    Population   = EnvironmentalSelection(Population,V,(Global.gen/Global.maxgen)^2);
    
    %% Optimization
    while Global.NotTermination(Population)
        Fitness = calFitness(Population.objs);
        if length(Population) >= 2
            Rank = randperm(length(Population),floor(length(Population)/2)*2);
        else
            Rank = [1,1];
        end
        Loser  = Rank(1:end/2);
        Winner = Rank(end/2+1:end);
        Change = Fitness(Loser) >= Fitness(Winner);
        Temp   = Winner(Change);
        Winner(Change) = Loser(Change);
        Loser(Change)  = Temp;
        Offspring      = Operator(Population(Loser),Population(Winner),Tasks);
        Population     = EnvironmentalSelection([Population,Offspring],V,(Global.gen/Global.maxgen)^2);
    end
end
function Offspring = Operator(Loser,Winner,Tasks)
% The competitive swarm optimizer of LMOCSO
    %% Parameter setting
    count=1;
    for i = 1:length(Winner)
        LoserDec  = Loser(i).dec;
        WinnerDec = Winner(i).dec;
        D1      = length(LoserDec);
        D2      = length(WinnerDec);
        Loser(i).adds2(zeros(1,D1));
        Winner(i).adds2(zeros(1,D2));
        if Loser(i).add == Winner(i).add
            D = D1;
            LoserVel =  Loser(i).add2;
            WinnerVel = Winner(i).add2;
        else
            D =  Tasks(Loser(i).add).D_high;
            LoserDec  = (Tasks(Loser(i).add).A*LoserDec')';
            WinnerDec = (Tasks(Winner(i).add).A*WinnerDec')';
            LoserVel =  (Tasks(Loser(i).add).A*Loser(i).add2')';
            WinnerVel = (Tasks(Winner(i).add).A*Winner(i).add2')';
        end
        % Competitive swarm optimizer
        r1 = rand(1,D);
        r2 = rand(1,D);
        OffVel = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
        OffDec = LoserDec + OffVel + r1.*(OffVel-LoserVel);
        % Polynomial mutation
        OffDec = mutate(OffDec,20,1);
        WinnerDec = mutate(WinnerDec,20,1);
        skill_factor = Winner(i).add;
        if D ~=  Tasks(skill_factor).D_eff
            OffDec = (Tasks(skill_factor).A_inv*OffDec')';
            OffVel = (Tasks(skill_factor).A_inv*OffVel')';
            WinnerDec = (Tasks(skill_factor).A_inv*WinnerDec')';
            WinnerVel = (Tasks(skill_factor).A_inv*WinnerVel')';
        end
        Offspring(count) = Individual(OffDec,Tasks,skill_factor,OffVel);
        Offspring(count+1) = Individual(WinnerDec,Tasks,skill_factor,WinnerVel);
        count=count+2;
    end
    

 
    
end
function Fitness = calFitness(PopObj)
% Calculate the fitness by shift-based density
    N      = size(PopObj,1);
    fmax   = max(PopObj,[],1);
    fmin   = min(PopObj,[],1);
    PopObj = (PopObj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
    Dis    = inf(N);
    for i = 1 : N
        SPopObj = max(PopObj,repmat(PopObj(i,:),N,1));
        for j = [1:i-1,i+1:N]
            Dis(i,j) = norm(PopObj(i,:)-SPopObj(j,:));
        end
    end
    Fitness = min(Dis,[],2);
end

function Population = EnvironmentalSelection(Population,V,theta)
% The environmental selection of LMOCSO
    Population = Population(NDSort(Population.objs,1)==1);
    PopObj = Population.objs;
    [N,M]  = size(PopObj);
    NV     = size(V,1);
    
    %% Translate the population
    PopObj = PopObj - repmat(min(PopObj,[],1),N,1);
    
    %% Calculate the smallest angle value between each vector and others
    cosine = 1 - pdist2(V,V,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    gamma  = min(acos(cosine),[],2);

    %% Associate each solution to a reference vector
    Angle = acos(1-pdist2(PopObj,V,'cosine'));
    [~,associate] = min(Angle,[],2);

    %% Select one solution for each reference vector
    Next = zeros(1,NV);
    for i = unique(associate)'
        current = find(associate==i);
        % Calculate the APD value of each solution
        APD = (1+M*theta*Angle(current,i)/gamma(i)).*sqrt(sum(PopObj(current,:).^2,2));
        % Select the one with the minimum APD value
        [~,best] = min(APD);
        Next(i) = current(best);
    end
    % Population for next generation
    Population = Population(Next(Next~=0));
end