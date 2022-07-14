function result = function2(x,y)
%exp r + 25 · sin (2.5 − r) · exp[− (2.5 − r)^2]
    result = [];
    for pt_x = x
        res = [];
        for pt_y = y
            r = sqrt(pt_x*pt_x + pt_y*pt_y);
            val = exp(r) + 25 * sin(2.5-r) * exp(-1*(2.5 - r)*(2.5 - r));
            res=[res,val];
        end
        result = [result;res];
    end
end

