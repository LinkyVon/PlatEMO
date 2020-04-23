function Population = subCMOPSO(Population,Tasks,rmp)
    Global  = GLOBAL.GetObj();
    while Global.NotTermination(Population)
        Offspring  = Operator(Population,Tasks);
        Population = EnvironmentalSelection([Population,Offspring],Global.N);
    end
end
function Offspring = Operator(Population,Tasks)
% The particle swarm optimization in CMOPSO
    %% Get leaders 
    Front     = NDSort(Population.objs,inf);    
    [~,rank]  = sortrows([Front',-CrowdingDistance(Population.objs,Front)']);
    LeaderSet = rank(1:10);
    count=1;
    for i = 1:length(Population)    
        % Parameter setting
        P_Dec  = Population(i).dec;     
        D      = length(P_Dec); 
        P_Obj  = Population(i).obj;
        Population(i).adds2(zeros(1,D)); 
        % Competition according to the angle
        winner = LeaderSet(randperm(length(LeaderSet),2));
        P_winner = Population(winner);
        c1     = dot(P_Obj,P_winner(1).obj)/(norm(P_Obj)*norm(P_winner(1).obj));
        angle1 = rad2deg(acos(c1));
        c2     = dot(P_Obj,P_winner(2).obj)/(norm(P_Obj)*norm(P_winner(2).obj));
        angle2 = rad2deg(acos(c2));
        mask   = (angle1 > angle2);
        winner = ~mask.*winner(1) + mask.*winner(2);
        % Learning
        if Population(i).add == Population(winner).add
            V =  Population(i).add2;
            winner_Dec = Population(winner).dec;
        else
            D =  Tasks(Population(i).add).D_high;
            P_Dec  = (Tasks(Population(i).add).A*Population(i).dec')';
            V =  (Tasks(Population(i).add).A*Population(i).add2')';
            winner_Dec = (Tasks(Population(winner).add).A*Population(winner).dec')';
        end
        r1 = rand(1,D);
        r2 = rand(1,D);
        Off_V = r1.*V + r2.*(winner_Dec-P_Dec);
        Off_P = P_Dec + Off_V;
        %Polynomial mutation
        Off_P = mutate(Off_P,20,1);
        skill_factor = Population(winner).add;
        if D ~=  Tasks(skill_factor).D_eff
            Off_P = (Tasks(skill_factor).A_inv*Off_P')';
            Off_V = (Tasks(skill_factor).A_inv*Off_V')';
        end
        Offspring(count) = Individual(Off_P,Tasks,skill_factor,Off_V);
        count=count+1;
    end
end

function Population = EnvironmentalSelection(Population,N)
% The environmental selection of CMOPSO
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,N);
    Next = false(1,length(FrontNo));
    Next(FrontNo<MaxFNo) = true;
    
    PopObj = Population.objs;
    fmax   = max(PopObj(FrontNo==1,:),[],1);
    fmin   = min(PopObj(FrontNo==1,:),[],1);
    PopObj = (PopObj-repmat(fmin,size(PopObj,1),1))./repmat(fmax-fmin,size(PopObj,1),1);

    %% Select the solutions in the last front
    Last = find(FrontNo==MaxFNo);
    del  = Truncation(PopObj(Last,:),length(Last)-N+sum(Next));
    Next(Last(~del)) = true;
    % Population for next generation
    Population = Population(Next);
end
function Del = Truncation(PopObj,K)
% Select part of the solutions by truncation

    N = size(PopObj,1);
    
    %% Truncation
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,N);
    while sum(Del) < K
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
end