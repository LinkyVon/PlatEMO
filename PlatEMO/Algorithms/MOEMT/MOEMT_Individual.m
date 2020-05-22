classdef MOEMT_Individual < handle
%INDIVIDUAL - The class of an individual.


    properties(SetAccess = public)
        dec;        % Decision variables of the individual
        inh;        % the actual individual with high dimension        
        obj;        % Objective values of the individual
        con;        % Constraint violations of the individual
        skill_factor;        % Additional properties of the individual
        rank;
        vel;        % Additional properties of the individual
    end
    methods
        %% Constructor
        function obj = MOEMT_Individual(Decs,Tasks,skill_factor,Vel)
        %INDIVIDUAL - Constructor of INDIVIDUAL class.
            if nargin > 0
                % Create new objects
                if length(Decs) == Tasks(skill_factor).D_func
                    Decs_high = (Tasks(skill_factor).A*Decs')';
                else
                    Decs_high = Decs;
                    Decs = (Tasks(skill_factor).A_inv*Decs')';
                end
                obj(1,size(Decs,1)) = MOEMT_Individual;
                Global = GLOBAL.GetObj();
                 
                % Set the infeasible decision variables to boundary values
                if ~isempty(Global.lower) && ~isempty(Global.upper)
                    Lower = repmat(Global.lower,length(obj),1);
                    Upper = repmat(Global.upper,length(obj),1);
                    Decs_high  = max(min(Decs_high,Upper),Lower);
                    Lower2 = repmat(Global.lower(1:length(Decs)),length(obj),1);
                    Upper2 = repmat(Global.upper(1:length(Decs)),length(obj),1);
                    Decs  = max(min(Decs,Upper2),Lower2);
                end
                % Calculte the objective values and constraint violations
                obj.dec = Decs;
                obj.inh = INDIVIDUAL(Decs_high);
                obj.obj = obj.inh.obj;
                obj.con = obj.inh.con;
                obj.skill_factor = skill_factor;
                
                if nargin > 3
                    obj.vel = Vel;
                end
                % Increase the number of evaluated individuals
                Global.evaluated = Global.evaluated + length(obj);
            end
        end
        %% Get the matrix of decision variables of the population
        function value = decs(obj)
        %decs - Get the matrix of decision variables of the population.
        %
        %   A = obj.decs returns the matrix of decision variables of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
            Decs = [];
            for i = 1:length(obj)
                Decs(i,:) = obj(i).inh.dec;
            end
            value = cat(1,Decs);
        end
        
        %% Get the matrix of objective values of the population
        function value = objs(obj)
        %objs - Get the matrix of objective values of the population.
        %
        %   A = obj.objs returns the matrix of objective values of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.obj);
        end
        %% Get the matrix of constraint violations of the population
        function value = cons(obj)
        %cons - Get the matrix of constraint violations of the population.
        %
        %   A = obj.cons returns the matrix of constraint violations of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
        
            value = cat(1,obj.con);
        end
        %% Get the matrix of additional properties of the population
        function value = skill_factors(obj)
        %adds - Get the matrix of additional properties of the population.
        %
        %   A = obj.adds(AddProper) returns the matrix of additional
        %   properties of the population obj. If any individual in obj does
        %   not contain an additional property, assign it a default value
        %   specified in AddProper.

            value = cat(1,obj.skill_factor);
        end
        %% Get vel
        function value = vels(obj,AddProper)
        %adds - Get the matrix of additional properties of the population.
        %
        %   A = obj.adds(AddProper) returns the matrix of additional
        %   properties of the population obj. If any individual in obj does
        %   not contain an additional property, assign it a default value
        %   specified in AddProper.

            for i = 1 : length(obj)
                if isempty(obj(i).vel)
                    obj(i).vel = AddProper(i,:);
                end
            end
            value = cat(1,obj.vel);
        end
        %% Get rank
        function value = ranks(obj,AddProper)
        %adds - Get the matrix of additional properties of the population.
        %
        %   A = obj.adds(AddProper) returns the matrix of additional
        %   properties of the population obj. If any individual in obj does
        %   not contain an additional property, assign it a default value
        %   specified in AddProper.

            for i = 1 : length(obj)
                if isempty(obj(i).rank)
                    obj(i).rank = AddProper(i,:);
                end
            end
            value = cat(1,obj.rank);
        end
    end
end