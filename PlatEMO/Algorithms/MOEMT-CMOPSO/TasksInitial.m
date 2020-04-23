function Tasks = TasksInitial(no_of_tasks,D,D_high)
    for i=1:no_of_tasks
        Tasks(i).D_eff = D(i);
        Tasks(i).D_high = D_high;
        if i ~= no_of_tasks
            Tasks(i).A = normrnd(0,1,Tasks(i).D_high,Tasks(i).D_eff);
        else
            Tasks(i).A = eye(D_high, D_high);
        end
        Tasks(i).A_inv = pinv(Tasks(i).A);  %(A'*A)\A'
        Tasks(i).Pop = 0;
    end
end