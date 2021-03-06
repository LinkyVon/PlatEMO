function Population = subLMOCSO(Population,Tasks,rmp)
    Global  = GLOBAL.GetObj();
%     [V,Global.N] = UniformPoint(Global.N,Global.M);
%     Population   = EnvironmentalSelection(Population,V,(Global.gen/Global.maxgen)^2);
    
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
        Offspring      = Operator(Population(Loser),Population(Winner),Tasks,rmp);
        Population = EnvironmentalSelection([Population,Offspring],Global.N);
%         Population     = EnvironmentalSelection([Population,Offspring],V,(Global.gen/Global.maxgen)^2);
        if mod(Global.gen,10) == 0 && length(Population)>Global.N/2
            Population = UpdateTask(Population,Tasks); 
        end
    end
end
function Offspring = Operator(Loser,Winner,Tasks,rmp)
% The competitive swarm optimizer of LMOCSO
    %% Parameter setting
    count=1;
    for i = 1:length(Winner)
        LoserDec  = Loser(i).dec;
        WinnerDec = Winner(i).dec;
        D1      = length(LoserDec);
        D2      = length(WinnerDec);
        Loser(i).vels(zeros(1,D1));
        Winner(i).vels(zeros(1,D2));
        if Loser(i).skill_factor == Winner(i).skill_factor
            D = D1;
            LoserVel =  Loser(i).vel;
            WinnerVel = Winner(i).vel;
        else
            D =  Tasks(Loser(i).skill_factor).D_high;
            LoserDec  = Loser(i).inh.dec; %(Tasks(Loser(i).skill_factor).A*LoserDec')';
            WinnerDec = Winner(i).inh.dec; %(Tasks(Winner(i).skill_factor).A*WinnerDec')';
            LoserVel =  (Tasks(Loser(i).skill_factor).A*Loser(i).vel')';
            WinnerVel = (Tasks(Winner(i).skill_factor).A*Winner(i).vel')';
        end
        % Competitive swarm optimizer
        r1 = rand(1,D);
        r2 = rand(1,D);
        if rand(1)<rmp || Loser(i).skill_factor == Winner(i).skill_factor
            OffVel = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
            if rand(1)<0.5
                skill_factor= Winner(i).skill_factor;
            else
                skill_factor= Loser(i).skill_factor;
            end
        else
            for j = 1: length(Winner)
                if Winner(j).skill_factor ~= Loser(i).skill_factor
                    WinnerDecDiff = Winner(j).inh.dec; %(Tasks(Winner(j).skill_factor).A*Winner(j).dec')';
                    r3 = rand(1,D);
                    OffVel = r1.*LoserVel + r2.*(WinnerDec-LoserDec)+r3.*(WinnerDecDiff-LoserDec);
                    if rand(1)<0.5
                        skill_factor= Winner(j).skill_factor;
                    else
                        skill_factor= Loser(i).skill_factor;
                    end
                    break;
                else
                    OffVel = r1.*LoserVel + r2.*(WinnerDec-LoserDec);
                    if rand(1)<0.5
                        skill_factor= Winner(i).skill_factor;
                    else
                        skill_factor= Loser(i).skill_factor;
                    end
                end
            end
        end
        OffDec = LoserDec + OffVel + r1.*(OffVel-LoserVel);
        % Polynomial mutation
        OffDec = mutate(OffDec,20,1);
        WinnerDec = mutate(WinnerDec,20,1);
        if D ~=  Tasks(skill_factor).D_func
            %OffDec = (Tasks(skill_factor).A_inv*OffDec')';
            OffVel = (Tasks(skill_factor).A_inv*OffVel')';
            %WinnerDec = (Tasks(skill_factor).A_inv*WinnerDec')';
            WinnerVel = (Tasks(skill_factor).A_inv*WinnerVel')';
        end
        Offspring(count) = MOEMT_Individual(OffDec,Tasks,skill_factor,OffVel);
        Offspring(count+1) = MOEMT_Individual(WinnerDec,Tasks,skill_factor,WinnerVel);
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

% function Population = EnvironmentalSelection(Population,V,theta)
% % The environmental selection of LMOCSO
%     Population = Population(NDSort(Population.objs,1)==1);
%     PopObj = Population.objs;
%     [N,M]  = size(PopObj);
%     NV     = size(V,1);
%     
%     %% Translate the population
%     PopObj = PopObj - repmat(min(PopObj,[],1),N,1);
%     
%     %% Calculate the smallest angle value between each vector and others
%     cosine = 1 - pdist2(V,V,'cosine');
%     cosine(logical(eye(length(cosine)))) = 0;
%     gamma  = min(acos(cosine),[],2);
% 
%     %% Associate each solution to a reference vector
%     Angle = acos(1-pdist2(PopObj,V,'cosine'));
%     [~,associate] = min(Angle,[],2);
% 
%     %% Select one solution for each reference vector
%     Next = zeros(1,NV);
%     for i = unique(associate)'
%         current = find(associate==i);
%         % Calculate the APD value of each solution
%         APD = (1+M*theta*Angle(current,i)/gamma(i)).*sqrt(sum(PopObj(current,:).^2,2));
%         % Select the one with the minimum APD value
%         [~,best] = min(APD);
%         Next(i) = current(best);
%     end
%     % Population for next generation
%     Population = Population(Next(Next~=0));
% end