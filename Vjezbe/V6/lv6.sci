function[x] = halfinterval(f, a, b, eps)
    f1 = f(a);
    f2 = f(b);
    pronadjeno = 0;
    while (b - a > eps)
        c = (a + b) / 2;
        f3 = f(c);
        if f3 == 0
            x = c;
            pronadjeno = 1;
            break;
        end
        if f1 * f3 <= 0
            b = c;
            f2 = f3;
        else
            a = c;
            f1 = f3;
        end           
    end
if pronadjeno ~= 1
    x = (a + b) / 2;
end     
endfunction;