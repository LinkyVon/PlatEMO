 function Population = Individual(Decs,Tasks,skill_factor,Vel)
    %INDIVIDUAL(obj.problem.Init(N));
    if nargin > 0
        % Create new objects
        Population(1,size(Decs,1)) = INDIVIDUAL;
        Global = GLOBAL.GetObj();
        % Set the infeasible decision variables to boundary values
        if Tasks(skill_factor).D_eff ~= Global.D
            Decs_high = (Tasks(skill_factor).A*Decs')';
        else
            Decs_high = Decs;
        end
        if ~isempty(Global.lower) && ~isempty(Global.upper)
            Lower = repmat(Global.lower,length(Population),1);
            Upper = repmat(Global.upper,length(Population),1);
            Decs_high  = max(min(Decs_high,Upper),Lower);
            Lower2 = repmat(Global.lower(1:length(Decs)),length(Population),1);
            Upper2 = repmat(Global.upper(1:length(Decs)),length(Population),1);
            Decs  = max(min(Decs,Upper2),Lower2);
        end
        % Calculte the objective values and constraint violations
        Decs = Global.problem.CalDec(Decs);
        Objs = Global.problem.CalObj(Decs_high);
        Cons = Global.problem.CalCon(Decs);
        % Assign the decision variables, objective values,
        % constraint violations, and additional properties
        for i = 1 : length(Population)
            Population(i).dec = Decs(i,:);
            Population(i).obj = Objs(i,:);
            Population(i).con = Cons(i,:);
        end
        if nargin > 2
            for i = 1 : length(Population)
                Population(i).add = skill_factor(i,:);
            end
            if nargin > 3
                for i = 1 : length(Population)
                    Population(i).add2 = Vel(i,:);
                end
            end
        end
        % Increase the number of evaluated individuals
        Global.evaluated = Global.evaluated + length(Population);
    end
end