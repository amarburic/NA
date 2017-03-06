function [s] = Lagrange(xdata1, ydata1, x)
    n1 = size(x, 2);
    s=zeros(1, n1);
    n2 = size(xdata1, 2);
    for i = 1:n2
        yi = ydata1(i);
        if yi == 0
            continue;
        end
        p = zeros(1, n1);
        p(1, 1:n1)=p(1, 1:n1) + yi;
        for j = 1:n2
            if(j<>i)
                p(1, 1:n1)=p(1, 1:n1). * (x(1, 1:n1) - xdata1(j))/(xdata1(i)-xdata1(j));
            end            
        end
        s=s+p;
    end
endfunction

function my_plot(x_min, x_max, n1, n2, fstring)
    xdata1 = linspace(x_min, x_max, n1);
    deff('[y]=f(x)','y = ' + fstring);
    ydata1 = feval(xdata1, f);
    x1 = linspace(x_min, x_max, n2);
    y1 = Lagrange(xdata1, ydata1, x1);
    plot(x1, y1);
endfunction

subplot(321)
my_plot(0, 2 * %pi, 10, 50, 'sin(x)');

subplot(322)
my_plot(-1, 1, 10, 50, '1 / (1 + 25 * x ^ 2)');

subplot(323)
my_plot(1, 10, 10, 50, 'log(x)');

subplot(324)
my_plot(-30, 30, 50, 100, 'x^3 * sin(x)');

//weierstrassova funkcija koja se racuna do n-tog sabirka
function [y] = w(x, a, b, n)
    y = 0;
    for i = 0:n
        y = y + a^i * cos(b^i * %pi * x);
    end
endfunction

subplot(325)
my_plot(-2, 2, 10, 50, 'w(x, 0.816, 7, 50)');



