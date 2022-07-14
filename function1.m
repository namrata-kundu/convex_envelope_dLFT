function result = function1(x,y)
% %(x^2 + y^2 - 1)^2

    result = [];
    for pt_x = x
        res = [];
        for pt_y = y
            val = (pt_x*pt_x + pt_y*pt_y - 1)*(pt_x*pt_x + pt_y*pt_y - 1);
            res=[res,val];
        end
        result = [result;res];
    end
end

