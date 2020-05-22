function Population = UpdateTask(Population,Tasks)
    % Update Task
    Pop = zeros(1,length(Tasks));
    for i = 1:length(Tasks)
        Pop(i) = sum([Population.skill_factor]==i);
    end
    for i = 1:length(Tasks)
        if Pop(i) <= floor(length(Population)/(2*length(Tasks))) && max(Pop)>= floor(length(Population)/3) 
            Tasks(i).A = normrnd(0,1,Tasks(i).D_high,Tasks(i).D_func);
            Tasks(i).A_inv = pinv(Tasks(i).A); 
            sf = find(Pop ==max(Pop)); 
            P = find(Population.skill_factors==sf(1));
            P = sort(P,'descend');
            for j = 1:floor(length(Population)/(2*length(Tasks)))
                Population(P(j)).skill_factor = i;
                Population(P(j)).dec = (Tasks(i).A_inv*Population(P(j)).inh.dec')';
                if ~isempty(Population(P(j)).vel)
                    Population(P(j)).vel = (Tasks(i).A_inv*(Tasks(sf(1)).A*Population(P(j)).vel'))';
                end
            end
        end
    end
end