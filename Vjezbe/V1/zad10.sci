function a(izraz)
    izraz = "y = " + izraz;
    deff("[y] = f(x)", izraz);
    x = linspace(-2, 2, 100);
    y = feval(x, f);
    subplot(121);
    plot(x, y, 'ob');
    subplot(122);
    plot(x, y, 'xr');
endfunction

function b(izraz)
    izraz = "z = " + izraz;
    deff("[z] = f(x, y)", izraz);
    x = linspace(-4, 4, 100);
    y = linspace(-4, 4, 100);
    z = feval(x, y, f);
    surf(x, y, z);
endfunction

