function result = function3(x,y)
% f1(x,y) + (for i=1 to 7 summation(f1(0,10i).e^(-10*(r-10i)*(r-10i)))

    result = [];
    for pt_x = x
        res = [];
        for pt_y = y
            val1 = function1(pt_x,pt_y);
            val2 = 0;
            for i=1:7
                r = sqrt(pt_x*pt_x + pt_y*pt_y);
                second_part = exp(-10*(r - 10*i)*(r - 10*i));
                val2 = val2 + (function1(0,10*i)*second_part);
            end
            res=[res,val1+val2];
        end
        result = [result;res];
    end
end

