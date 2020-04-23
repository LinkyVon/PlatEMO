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
        if parent(p1).add==parent(p2).add
            [Decs_1,Decs_2] = crossover(parent(p1).dec,parent(p2).dec,disC,proC); 
            Decs_1 = mutate(Decs_1,disM,proM);
            Decs_2 = mutate(Decs_2,disM,proM);
            if rand(1)<0.5
                skill_factor= [parent(p1).add,parent(p1).add];
            else
                skill_factor= [parent(p2).add,parent(p2).add];
            end
        else
            if rand(1)<rmp % different skill_factor
                P1 = (Tasks(parent(p1).add).A*parent(p1).dec')';
                P2 = (Tasks(parent(p2).add).A*parent(p2).dec')';
                [Decs_1,Decs_2] = crossover(P1,P2,disC,proC); 
                Decs_1 = mutate(Decs_1,disM,proM);
                Decs_2 = mutate(Decs_2,disM,proM);
                if rand(1)<0.5
                    skill_factor= [parent(p1).add,parent(p1).add];
                    Decs_1 = (Tasks(parent(p1).add).A_inv*Decs_1')';
                    Decs_2 = (Tasks(parent(p1).add).A_inv*Decs_2')';
                else
                    skill_factor= [parent(p2).add,parent(p2).add];
                    Decs_1 = (Tasks(parent(p2).add).A_inv*Decs_1')';
                    Decs_2 = (Tasks(parent(p2).add).A_inv*Decs_2')';
                end
            else
                Decs_1 = mutate(parent(p1).dec,disM,proM);
                Decs_2 = mutate(parent(p2).dec,disM,proM);
                skill_factor= [parent(p1).add,parent(p2).add];
            end
        end
        Offspring(count) = Individual(Decs_1,Tasks,skill_factor(1));
        Offspring(count+1) = Individual(Decs_2,Tasks,skill_factor(2));
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
% function Offspring = SBX(Parent,Parameter)
%     %% Parameter setting
%     if nargin > 1
%         [proC,disC] = deal(Parameter{:});
%     else
%         [proC,disC] = deal(1,20);
%     end
%     [N,D]   = size(Parent(1));
%     beta = zeros(N,D);
%     mu   = rand(N,D);
%     beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
%     beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
%     beta = beta.*(-1).^randi([0,1],N,D);
%     beta(rand(N,D)<0.5) = 1;
%     beta(repmat(rand(N,1)>proC,1,D)) = 1;
%     Offspring = [(Parent(1)+Parent(2))/2+beta.*(Parent(1)-Parent(2))/2
%                  (Parent(1)+Parent(2))/2-beta.*(Parent(1)-Parent(2))/2];
% end
% function Offspring = mutate(Decs,Parameter)
%     %% Parameter setting
%     if nargin > 1
%         [proM,disM] = deal(Parameter{:});
%     else
%         [proM,disM] = deal(1,20);
%     end
%     [N,D]   = size(Decs);
%     Global  = GLOBAL.GetObj();
%     Lower = repmat(Global.lower(1:D),1*N,1);
%     Upper = repmat(Global.upper(1:D),1*N,1);
%     Site  = rand(1*N,D) < proM/D;
%     mu    = rand(1*N,D);
%     temp  = Site & mu<=0.5;
%     Offspring       = min(max(Decs,Lower),Upper);
%     Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%                       (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
%     temp = Site & mu>0.5; 
%     Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%                       (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
% end