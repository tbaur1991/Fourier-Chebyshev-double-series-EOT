function Ts = chebyshevList(nu,u)
    if nu > 0
        % first 10 Chebyshev polynomials
        Ts = zeros(1,10);
        if ~isnumeric(u)
            Ts = subs(Ts);
        end
        Ts(1) = u;
        Ts(2) = 2*u^2 - 1;
        Ts(3) = 4*u^3 - 3*u;
        Ts(4) = 8*u^4 - 8*u^2 + 1;
        Ts(5) = 16*u^5 - 20*u^3 + 5*u;
        Ts(6) = 32*u^6 - 48*u^4 + 18*u^2 - 1;
        Ts(7) = 64*u^7 - 112*u^5 + 56*u^3 - 7*u;
        Ts(8) = 128*u^8 - 256*u^6 + 160*u^4 - 32*u^2 + 1;
        Ts(9) = 256*u^9 - 576*u^7 + 432*u^5 - 120*u^3 + 9*u;
        Ts(10) = 512*u^10 - 1280*u^8 + 1120*u^6 - 400*u^4 + 50*u^2 - 1;
        % output 10 Chebyshev polynomials
        Ts = Ts(1:nu);
    else
        Ts = 1;
    end
end