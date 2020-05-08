classdef INDIVIDUAL1 < handle
%INDIVIDUAL - The class of an individual.
%
%   This is the class of an individual. An object of INDIVIDUAL stores all
%   the properties including decision variables, objective values,
%   constraint violations, and additional properties of an individual.
%
% INDIVIDUAL properties:
%   dec         <public>     decision variables of the individual
%   dec_high         <public>     high decision variables of the individual
%   obj         <public>     objective values of the individual
%   con         <public>     constraint violations of the individual
%   skill_factor         <public>     additional properties of the individual
%   vel         <public>     additional properties of the individual
%
% INDIVIDUAL methods:
%   INDIVIDUAL	<public>        the constructor, all the properties will be
%                               set when the object is creating
%   decs        <public>      	get the matrix of high decision variables of the
%                               population
%   dec_highs        <public>      	get the matrix of decision variables of the
%                               population
%   objs        <public>        get the matrix of objective values of the
%                               population
%   cons        <public>        get the matrix of constraint violations of
%                               the population
%   skill_factors        <public>        get the matrix of additional properties of
%                               the population
%   vel        <public>       get the matrix of additional properties of
%                               the population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = public)
        dec;        % Decision variables of the individual
        dec_high;        % high Decision variables of the individual        
        obj;        % Objective values of the individual
        con;        % Constraint violations of the individual
        skill_factor;        % Additional properties of the individual
        vel;        % Additional properties of the individual
    end
    methods
        %% Constructor
        function obj = INDIVIDUAL1(Decs,Tasks,skill_factor,Dec_high,Vel)
        %INDIVIDUAL - Constructor of INDIVIDUAL class.
        %
        %   H = INDIVIDUAL(Dec) creates an array of individuals (i.e., a
        %   population), where Dec is the matrix of decision variables of
        %   the population. The objective values and constraint violations
        %   are automatically calculated by the test problem functions.
        %   After creating the individuals, the number of evaluations will
        %   be increased by length(H).
        %
        %   H = INDIVIDUAL(Dec,AddProper) creates the population with
        %   additional properties stored in AddProper, such as the velocity
        %   in particle swarm optimization.
        %
        %   Example:
        %       H = INDIVIDUAL(rand(100,3))
        %       H = INDIVIDUAL(rand(100,10),randn(100,3))
        
            if nargin > 0
                % Create new objects
                obj(1,size(Decs,1)) = INDIVIDUAL1;
                Global = GLOBAL.GetObj();
                if skill_factor == 1
                    Decs_high = Dec_high;
                elseif skill_factor == length(Tasks)
                    Decs_high = Decs;
                else
                    Decs_high = (Tasks(skill_factor).A*Decs')';
                end
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
                
                Z = transformation(Decs,Decs_high,skill_factor);
                Objs = Global.problem.CalObj(Z);
                Decs = Global.problem.CalDec(Decs);
                Cons = Global.problem.CalCon(Decs);
                % Assign the decision variables, objective values,
                % constraint violations, and additional properties
                for i = 1 : length(obj)
                    obj(i).dec = Decs(i,:);
                    obj(i).obj = Objs(i,:);
                    obj(i).con = Cons(i,:);
                    obj(i).skill_factor = skill_factor(i,:);
                    obj(i).dec_high = Decs_high(i,:);
                end
                if nargin > 4
                    for i = 1 : length(obj)
                       obj(i).vel = Vel(i,:);
                    end
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
        %change
            Global = GLOBAL.GetObj();
            Decs = [];
            for i = 1:length(obj)
                if length(obj(i).dec) ~= Global.D
                    t = [obj(i).dec, zeros(1,Global.D-length(obj(i).dec))];
                    Decs(i,:) = t;
                else
                    Decs(i,:) = obj(i).dec;
                end
            end

            value = cat(1,Decs);
        end
         %% Get the matrix of decision variables of the population
        function value = dec_highs(obj)
        %decs - Get the matrix of decision variables of the population.
        %
        %   A = obj.decs returns the matrix of decision variables of the
        %   population obj, where obj is an array of INDIVIDUAL objects.
       
            value = cat(1,obj.dec_high);
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
        function value = skill_factors(obj,AddProper)
        %adds - Get the matrix of additional properties of the population.
        %
        %   A = obj.adds(AddProper) returns the matrix of additional
        %   properties of the population obj. If any individual in obj does
        %   not contain an additional property, assign it a default value
        %   specified in AddProper.

            for i = 1 : length(obj)
                if isempty(obj(i).skill_factor)
                    obj(i).skill_factor = AddProper(i,:);
                end
            end
            value = cat(1,obj.skill_factor);
        end
        
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
    end
end