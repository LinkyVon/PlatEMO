function Offspring = OffspringCreation(parent,Tasks,rmp,Parameter)
    %% Parameter setting
    if nargin > 3
        [proC,disC,proM,disM] = deal(Parameter{:});
    else
        [proC,disC,proM,disM] = deal(1,20,1,20);
%        [proC,disC,proM,disM] = deal(0,10,1,10);
    end

    count=1;
    for i=1:2:length(parent)
        p1=i;
        p2=i+1;
        if parent(p1).skill_factor==parent(p2).skill_factor
            [Decs_1,Decs_2] = crossover(parent(p1).dec,parent(p2).dec,disC,proC); 
            Decs_1 = mutate(Decs_1,disM,proM);
            Decs_2 = mutate(Decs_2,disM,proM);
            if rand(1)<0.5
                skill_factor= [parent(p1).skill_factor,parent(p1).skill_factor];
            else
                skill_factor= [parent(p2).skill_factor,parent(p2).skill_factor];
            end
        else
            if rand(1)<rmp % different skill_factor
                P1 = parent(p1).inh.dec; %(Tasks(parent(p1).skill_factor).A*parent(p1).dec')';
                P2 = parent(p2).inh.dec; %(Tasks(parent(p2).skill_factor).A*parent(p2).dec')';
                [Decs_1,Decs_2] = crossover(P1,P2,disC,proC); 
                Decs_1 = mutate(Decs_1,disM,proM);
                Decs_2 = mutate(Decs_2,disM,proM);
                if rand(1)<0.5
                    skill_factor= [parent(p1).skill_factor,parent(p1).skill_factor];
                else
                    skill_factor= [parent(p2).skill_factor,parent(p2).skill_factor];
                end
            else
                Decs_1 = mutate(parent(p1).dec,disM,proM);
                Decs_2 = mutate(parent(p2).dec,disM,proM);
                skill_factor= [parent(p1).skill_factor,parent(p2).skill_factor];
            end
        end
        Offspring(count) = MOEMT_Individual(Decs_1,Tasks,skill_factor(1));
        Offspring(count+1) = MOEMT_Individual(Decs_2,Tasks,skill_factor(2));
        count=count+2;
    end     
end
function [rnvec1,rnvec2]=crossover(pvec1,pvec2,muc,prob_vswap)
    dim = length(pvec1);
    rnvec1=zeros(1,dim);
    rnvec2=zeros(1,dim);
    randlist=rand(1,dim);
    for i=1:dim
        if randlist(i)<=0.5
            k=(2*randlist(i))^(1/(muc+1));
        else
            k=(1/(2*(1-randlist(i))))^(1/(muc+1));
        end
        rnvec1(i)=0.5*(((1 + k)*pvec1(i)) + (1 - k)*pvec2(i));
        rnvec2(i)=0.5*(((1 - k)*pvec1(i)) + (1 + k)*pvec2(i));
        if rand(1) < prob_vswap
            tmp = rnvec1(i);
            rnvec1(i) = rnvec2(i);
            rnvec2(i) = tmp;
        end
        if rnvec1(i)<0
            rnvec1(i)=0;
        elseif rnvec1(i)>1
            rnvec1(i)=1;
        end
        if rnvec2(i)<0
            rnvec2(i)=0;
        elseif rnvec2(i)>1
            rnvec2(i)=1;
        end
    end            
end
function rnvec=mutate(p,mum,prob_mut)
    rnvec=p;
    dim = length(p);
    for i=1:dim
        if rand(1)<prob_mut/dim
            u=rand(1);
            if u <= 0.5
                del=(2*u)^(1/(1+mum)) - 1;
                rnvec(i)=p(i) + del*(p(i));
            else
                del= 1 - (2*(1-u))^(1/(1+mum));
                rnvec(i)=p(i) + del*(1-p(i));
            end
        end
    end            
end
