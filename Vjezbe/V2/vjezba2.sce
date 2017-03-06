function[x] = lijevo_dijeljenje(a, b)
    [n, m] = size(b);    
    
    for k = 1:n
        p = k
        for i = (k + 1):n
            u = a(i, k) / a(k, k);
            a(i, (k + 1):n) = a(i, (k + 1):n) - u * a(k, (k + 1):n);
            b(i, 1:m) = b(i, 1:m) - u * b(k, 1:m);           
        end
    end
    x = zeros(n, m);
    for k = 1:m
        for i = n : -1 : 1
            s = b(i, k);
            s = s - a(i, (i + 1):n) * x((i + 1):n, k);
            x(i, k) = s / a(i, i); 
        end
        
    end
endfunction
