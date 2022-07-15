function result = maxval_except(array, except)
    result = max(max(array(1:except-1)), max(array(except+1:end)));
end

