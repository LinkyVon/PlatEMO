function Population = subCMOPSO(Population,Tasks,rmp)
    Global  = GLOBAL.GetObj();
    while Global.NotTermination(Population)
        Offspring  = Operator(Population,Tasks,rmp);
        Population = EnvironmentalSelection([Population,Offspring],Global.N);
        if mod(Global.gen,10) == 0 
            Population = UpdateTask(Population,Tasks); 
        end 
    end
end
function Offspring = Operator(Population,Tasks,rmp)
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
        Population(i).vels(zeros(1,D)); 
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
        if Population(i).skill_factor == Population(winner).skill_factor
            V =  Population(i).vel;
            winner_Dec = Population(winner).dec;
        else
            D =  Tasks(Population(i).skill_factor).D_high;
            P_Dec  = Population(i).inh.dec; %(Tasks(Population(i).skill_factor).A*Population(i).dec')';
            V =  (Tasks(Population(i).skill_factor).A*Population(i).vel')';
            winner_Dec = Population(winner).inh.dec; %(Tasks(Population(winner).skill_factor).A*Population(winner).dec')';
        end
        r1 = rand(1,D);
        r2 = rand(1,D);
        if rand(1)<rmp || Population(i).skill_factor == Population(winner).skill_factor
            Off_V = r1.*V + r2.*(winner_Dec-P_Dec);
            if rand(1)<0.5
                skill_factor= Population(winner).skill_factor;
            else
                skill_factor= Population(i).skill_factor;
            end
        else
            for j = 1: length(LeaderSet)
                leader = LeaderSet(j);
                if Population(leader).skill_factor ~= Population(i).skill_factor
                    leader_Dec = Population(leader).inh.dec; %(Tasks(Population(leader).skill_factor).A*Population(leader).dec')';
                    r3 = rand(1,D);
                    Off_V = r1.*V + r2.*(winner_Dec-P_Dec)+r3.*(leader_Dec-P_Dec);
                    if rand(1)<0.5
                        skill_factor= Population(leader).skill_factor;
                    else
                        skill_factor= Population(i).skill_factor;
                    end
                    break;
                else
                    Off_V = r1.*V + r2.*(winner_Dec-P_Dec);
                    if rand(1)<0.5
                        skill_factor= Population(winner).skill_factor;
                    else
                        skill_factor= Population(i).skill_factor;
                    end
                    break;
                end
            end
        end
        Off_P = P_Dec + Off_V;% + r1.*(Off_V-V);
        %Polynomial mutation
        Off_P = mutate(Off_P,20,1);
        if D ~=  Tasks(skill_factor).D_func
            %Off_P = (Tasks(skill_factor).A_inv*Off_P')';
            Off_V = (Tasks(skill_factor).A_inv*Off_V')';
        end
        Offspring(count) = MOEMT_Individual(Off_P,Tasks,skill_factor,Off_V);
        count=count+1;
    end
end
