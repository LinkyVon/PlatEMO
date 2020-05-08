function value = transformation(Decs,Decs_high,skill_factor)
    if skill_factor ~=1
        value = Decs_high;
    else
        Global = GLOBAL.GetObj();
        minVal = Global.lower;
        maxVal = Global.upper;
        varsPerGroup = floor(length(Decs_high)/length(Decs));
        weight = [];
        for i = 1:length(Decs)-1
           weight = [weight, ones(1,varsPerGroup).*Decs(i)];
        end
        weight = [weight, ones(1,length(Decs_high)-size(weight,2)).*Decs(end)];
        pWert = 0.2;
        value = Decs_high+pWert.*(weight-1.0).*(maxVal-minVal);
         %do repair
        for i = 1:length(value)
            if value(i) < minVal(i)
               value(i) = minVal(i);
            elseif value(i) > maxVal(i)
               value(i) = maxVal(i);
            end
        end
    end
end