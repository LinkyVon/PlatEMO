function Population = subSMPSO(Population,Tasks,rmp) 
    Global  = GLOBAL.GetObj();
    Pbest            = Population;
    [Gbest,CrowdDis] = UpdateGbest(Population,Global.N);
    
    %% Optimization
    while Global.NotTermination(Gbest)
        ParentPool = TournamentSelection(2,Global.N,-CrowdDis);
        Population = Operator(Population,Pbest,Gbest(ParentPool),Global,Tasks,rmp);
        [Gbest,CrowdDis] = UpdateGbest([Gbest,Population],Global.N);
        Pbest            = UpdatePbest(Pbest,Population);
        
    end
    
end

function Offspring = Operator(Particle,Pbest,Gbest,Global,Tasks,rmp)
% Particle swarm optimization in SMPSO
    count=1;
    for i = 1:length(Particle)
        Particle(i).vels(zeros(1,length(Particle(i).dec)));
        if Particle(i).skill_factor == Pbest(i).skill_factor && Particle(i).skill_factor == Gbest(i).skill_factor
            D =  Tasks(Particle(i).skill_factor).D_eff;
            ParticleDec = Particle(i).dec;
            PbestDec    = Pbest(i).dec;
            GbestDec    = Gbest(i).dec;
            ParticleVel = Particle(i).vel;
        else
            D =  Tasks(Particle(i).skill_factor).D_high;
            ParticleDec = transformation(Particle(i).dec,Particle(i).dec_high,Particle(i).skill_factor);
            PbestDec    = transformation(Pbest(i).dec,Pbest(i).dec_high,Pbest(i).skill_factor);
            GbestDec    = transformation(Gbest(i).dec,Gbest(i).dec_high,Gbest(i).skill_factor);
            ParticleVel = (Tasks(Particle(i).skill_factor).A*Particle(i).vel')';
        end
        
        %% Particle swarm optimization
        W  = repmat(unifrnd(0.1,0.5),1,D);
        r1 = repmat(rand(1),1,D);
        r2 = repmat(rand(1),1,D);
        r3 = repmat(rand(1),1,D);
        C1 = repmat(unifrnd(1.5,2.5),1,D);
        C2 = repmat(unifrnd(1.5,2.5),1,D);
        C3 = repmat(unifrnd(1.5,2.5),1,D);
        if rand(1)<rmp || Particle(i).skill_factor == Gbest(i).skill_factor
            OffVel = W.*ParticleVel + C1.*r1.*(PbestDec-ParticleDec) + C2.*r2.*(GbestDec-ParticleDec);
            phi    = max(4,C1+C2);
            skill_factor= Gbest(i).skill_factor;
        else
            flag = 0;
            for j = 1: length(Gbest)
                if Gbest(j).skill_factor ~= Particle(i).skill_factor
                    flag = 1;
                    GbestDecDiff = transformation(Gbest(j).dec,Gbest(j).dec_high,Gbest(j).skill_factor);
                    OffVel = W.*ParticleVel + C1.*r1.*(PbestDec-ParticleDec) + C2.*r2.*(GbestDec-ParticleDec)+ C3.*r3.*(GbestDecDiff-ParticleDec);
                    phi    = max(6,C1+C2+C3);
                    if rand(1)<0.5
                        skill_factor= Gbest(j).skill_factor;
                    else
                        skill_factor= Particle(i).skill_factor;
                    end
                    break;
                end
            end
            if flag == 0
                OffVel = W.*ParticleVel + C1.*r1.*(PbestDec-ParticleDec) + C2.*r2.*(GbestDec-ParticleDec);
                phi    = max(4,C1+C2);
                if rand(1)<0.5
                    skill_factor= Gbest(i).skill_factor;
                else
                    skill_factor= Particle(i).skill_factor;
                end
            end
            
        end
        OffVel = OffVel.*2./abs(2-phi-sqrt(phi.^2-4*phi));
        delta  = (Global.upper(1:D)-Global.lower(1:D))/2;
        OffVel = max(min(OffVel,delta),-delta);
        OffDec = ParticleDec + OffVel;

        %% Polynomial mutation
        OffDec = mutate(OffDec,20,1);
        if D ~=  Tasks(skill_factor).D_eff
            OffDec_eff = (Tasks(skill_factor).A_inv*OffDec')';
            OffVel = (Tasks(skill_factor).A_inv*OffVel')';
        elseif skill_factor == 1
            OffDec_eff = OffDec;
            OffDec = Gbest(i).dec_high;
        else
            OffDec_eff = OffDec;
        end
        Offspring(count) = INDIVIDUAL1(OffDec_eff,Tasks,skill_factor,OffDec,OffVel);
        count=count+1;
    end
end

function [Gbest,CrowdDis] = UpdateGbest(Gbest,N)
% Update the global best set
    Gbest    = Gbest(NDSort(Gbest.objs,1)==1);
    CrowdDis = CrowdingDistance(Gbest.objs);
    [~,rank] = sort(CrowdDis,'descend');
    Gbest    = Gbest(rank(1:min(N,length(Gbest))));
    CrowdDis = CrowdDis(rank(1:min(N,length(Gbest))));
end

function CrowdDis = CrowdingDistance(PopObj)
% Calculate the crowding distance of each solution in the same front
    [N,M]    = size(PopObj);  
    CrowdDis = zeros(1,N);
    Fmax     = max(PopObj,[],1);
    Fmin     = min(PopObj,[],1);
    for i = 1 : M
        [~,rank] = sortrows(PopObj(:,i));
        CrowdDis(rank(1))   = inf;
        CrowdDis(rank(end)) = inf;
        for j = 2 : N-1
            CrowdDis(rank(j)) = CrowdDis(rank(j))+(PopObj(rank(j+1),i)-PopObj(rank(j-1),i))/(Fmax(i)-Fmin(i));
        end
    end
end
function Pbest = UpdatePbest(Pbest,Population)
% Update the local best position of each particle
    replace        = ~all(Population.objs>=Pbest.objs,2);
    Pbest(replace) = Population(replace);
end