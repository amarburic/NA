function [u] = polyfit(x,y,n)
    m = size(x,2);
    A=ones(m,n+1);
    y=y.';
    for i = 1:m
        for j = 2:n+1
            A(i,j) = A(i,j-1) * x(i);
        end
    end
    u = pinv(A)*y;
endfunction

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

function [u] = polyfit2(x,y,n)
    m = size(x,2);
    A=ones(m,n+1);
    y=y.';
    for i = 1:m
        for j = 2:n+1
            A(i,j) = A(i,j-1) * x(i);
        end
    end
    u = lijevo_dijeljenje(A.' * A, A.' * y);
endfunction

function [u] = polyfit3(x,y,n)
    m = size(x,2);
    A=ones(m,n+1);
    y=y.';
    for i = 1:m
        for j = 2:n+1
            A(i,j) = A(i,j-1) * x(i);
        end
    end
    [Q, R, E] = qr(A);
    u = R \ (Q.' * y);
    u = E * u;
endfunction

function [y] = polyval(u,x,n)
    m = size(x,2);
    y = zeros(1, m);
    for i = 1:n+1
        y(1, 1:m) = y(1, 1:m) + u(i) * x(1, 1:m) .^(i-1)
    end
endfunction

function test(xmin, xmax, n_tacki, func, n_poli)    
    x = linspace(xmin, xmax, n_tacki);
    deff('[y]=f(x)', 'y=' + func);
    y = feval(x, f);
    u = polyfit3(x, y, n_poli);
    y2 = polyval(u, x, n_poli);
    plot(x, y);
    plot(x, y2, 'rx');
    l2 = 'Polinom ' + string(n_poli) + '. reda';
    legend(func, l2);
endfunction

subplot(311)
test(-10, 10, 100, 'sin(x)', 10)

subplot(312)
test(1, 5, 20, 'log2(x)', 5)

subplot(313)
test(-5, 5, 100, 'sin(1/x)', 10)
