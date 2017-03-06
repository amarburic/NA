function [I] = simpson(f,a,b,eps)
    Nmax = 10;
    while(1)
    Nmin = 2;
    N = 2;
    h = (b-2)/N;
    s = (f(a) + f(b))/2;
    Told = s;
    Iold = 0;
    while(N<Nmax)
        for(i = 1:1N/2))
            s = s+f(a+(2*i-1)*h);
        end
        I = (4*h*s-Told)/3;
        if((N>Nmin)*(abs(I-Iold)<eps))
            return;
        end
        Told = h*s;
        Iold = I;
        N = 2*N;
        h = h/2;
    end

function [I] = simpson(f,a,b,eps)
    Nmax = 10;
    while(1)
    Nmin = 2;
    N = 2;
    h = (b-a)/N;
    s = (f(a) + f(b))/2;
    Told = s;
    Iold = 0;
    while(N<Nmax)
        for(i = 1:1N/2))
            s = s+f(a+(2*i-1)*h);
        end
        I = (4*h*s-Told)/3;
        if((N>Nmin)*(abs(I-Iold)<=eps))
            return;
        end
        Told = h*s;
        Iold = I;
        N = 2*N;
        h = h/2;
    end    
    Nmax = Nmax *2;
    end
    
endfunction
