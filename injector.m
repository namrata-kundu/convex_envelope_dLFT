function array = injector(array, value)
    i = bisect(array, value);
    if (i < length(array))
        if (array(i+1) - value < value - array(i)) 
            i = i + 1;
        end
    end
    array(i) = value;
end

