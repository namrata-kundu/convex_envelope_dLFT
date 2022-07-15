function result = function4(x,y)
% % 0.5(cos^2theta +

    result = [];
    alpha = -0.5;
    theta = 0.5*atan(1/3);
    for pt_x = x
        res = [];
        for pt_y = y
            part_1 = 0.5*(cos(theta).^2 + alpha*(sin(theta).^2))*pt_x*pt_x;
            part_2 = (1-alpha)*cos(theta)*sin(theta*pt_x*pt_y);
            part_3 = 0.5*(alpha*(cos(theta).^2) + sin(theta).^2)*pt_y*pt_y;
            val = part_1 + part_2 * part_3;
            res=[res,val];
        end
        result = [result;res];
    end
end

