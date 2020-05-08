 function MOEMT_WOF(Global)
% <algorithm> <EMT>
% Large-scale multi-objective  optimization via evolutionary multitasking
% D_func ---  30 --- low effective dim
% no_of_tasks ---  5 --- Number of tasks
% rmp --- 0.8---random mating probability
% solver ---  3 --- 1.NSGA II, 2.MOEAD, 3.SMPSO, 4.CSO
    %% Parameter settings
    [D_func,no_of_tasks,rmp,solver] = Global.ParameterSet(30,5,1,1);
    
    %% Generate random population and divide it into tasks based on skill factor
    while mod(Global.N,no_of_tasks) ~= 0
        Global.N = Global.N + 1;
    end
    N = Global.N;
    D_high = Global.D;
    D=zeros(1,no_of_tasks);
    for i=1:no_of_tasks-1
        D(i)= D_func;
    end
    D(no_of_tasks)= D_high; %D = [D_func, D_func,..,D_high]
    Tasks = TasksInitial(no_of_tasks,D,D_high);
    skill_factor = 1;
    for i=1:N
        D_T = D(skill_factor);
        PopDec = unifrnd(Global.lower(1:D_T),Global.upper(1:D_T)); % problem.Init(N)
        if skill_factor == 1
            Dec_high = unifrnd(Global.lower(1:D_high),Global.upper(1:D_high));
            Population(i) = INDIVIDUAL1(PopDec,Tasks,skill_factor,Dec_high);
        else
            Population(i) = INDIVIDUAL1(PopDec,Tasks,skill_factor);
        end
        if mod(i,N/no_of_tasks) == 0
            skill_factor= skill_factor +1;
        end
    end
    
%     for i=1:no_of_tasks
%         P = find([Population.add]==i);
%         Tasks(i).Pop = Population(P);
%     end
    
    
    %% Optimization
    switch solver
        case 1
            Population = subNSGAII(Population,Tasks,rmp);
        case 2
            Population = subMOEAD(Population,Tasks,rmp);
        case 3
            Population = subSMPSO(Population,Tasks,rmp);
        case 4
            Population = subCMOPSO(Population,Tasks,rmp);
        case 5 
            Population = subLMOCSO(Population,Tasks,rmp);
    end
end




