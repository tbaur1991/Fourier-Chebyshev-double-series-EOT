function Ts = chebyshevListDeriv(nu,u)
    if nu > 0
        % first 10 Chebyshev polynomials derivatives
        Ts = zeros(1,10);
        Ts(1) = 1;
        Ts(2) = 4*u;
        Ts(3) = 12*u^2 - 3;
        Ts(4) = 32*u^3 - 16*u;
        Ts(5) = 80*u^4 - 60*u^2 + 5;
        Ts(6) = 192*u^5 - 192*u^3 + 36*u;
        Ts(7) = 448*u^6 - 560*u^4 + 168*u^2 - 7;
        Ts(8) = 1024*u^7 - 1536*u^5 + 640*u^3 - 64*u;
        Ts(9) = 2304*u^8 - 4032*u^6 + 2160*u^4 - 360*u^2 + 9;
        Ts(10) = 5120*u^9 - 10240*u^7 + 6720*u^5 - 1600*u^3 + 100*u;
        % output 10 Chebyshev polynomials derivatives
        Ts = Ts(1:nu);
    else
        Ts = 1;
    end
end