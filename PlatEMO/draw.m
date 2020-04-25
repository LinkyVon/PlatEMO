pf = PF1(3,100);
score = zeros(1,100);
for i = 1:100
   score(i) = IGD(result{i,2}.objs,pf);  
end
 
function P = PF1(M,N)
    P = UniformPoint(N,M);
end
function P = PF9(M,N)
    interval     = [0,0.251412,0.631627,0.859401];
    median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
    X            = ReplicatePoint(N,M-1);
    X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
    X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
    P            = [X,2*(M-sum(X/2.*(1+sin(3*pi.*X)),2))];
end
function W = ReplicatePoint(SampleNum,M)
    if M > 1
        SampleNum = (ceil(SampleNum^(1/M)))^M;
        Gap       = 0:1/(SampleNum^(1/M)-1):1;
        eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
        eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
    else
        W = (0:1/(SampleNum-1):1)';
    end
end
function Score = IGD(PopObj,PF)
    Distance = min(pdist2(PF,PopObj),[],2);
    Score    = mean(Distance);
end