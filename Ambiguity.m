function y = Ambiguity(x, t, f)

    P = length(t);  F = length(f);  % T,F
    dt = t(2) - t(1);   % dt
    df = f(2) - f(1);   % df
    
    mo = round(f/df);  % mo
%     no = round(t/dt);  % no

    N = round(1/dt/df);
    y = zeros(P, F); 

    index = find(x);
    n1 = index(1);
    n2 = index(end);
    
    
    for p = 1:P
        
        % in this case, 1001 is the place t = 0
        pp = p;
        if(p > 1001)
            pp = pp - 1001;
        else
            pp = 1001 - pp;
        end

        Q = n2 - n1 - 2 * pp;
        c1 = zeros(1, N);
        
        for n = 0 : Q
            if(p > 1001)
                c1(n + 1) = x(n + n1 + 2 * pp) * conj(x(n + n1));
            else
                c1(n + 1) = x(n + n1) * conj(x(n + n1 + 2 * pp));
            end
        end

        C = fft(c1);

        for m = 1:F
            y(p, m) = C(mod(mo(m), N)+1)*exp(-1i * 2 * pi * Q * mo(m)/N)*dt;
        end
    end
end