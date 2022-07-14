function index = bisect(array, x)
    if (length(array) == 1)
        index = 1;
        return
    end

    if (x == array(length(array))) 
        index = length(array);
        return
    end

    a = 1;
    b = length(array);
    while true
        index = floor((a + b) / 2);
        if (array(index) <= x && array(index + 1) > x)
            return;
        elseif (array(index) > x) 
            b = index;
        else
            a = index;
        end
    end
end

