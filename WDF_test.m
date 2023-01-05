% ==================================================
% Write WDF to help check the result for CohenClass
% ==================================================
% The two function s and r is from p.96
% ==================================================
t = -9:0.01:1;
s = exp(1i * t.^2 / 10 - 1i * 3 * t);
s = [zeros(1, 100) s zeros(1, 900)];

t = -10:0.01:10;
r = exp(1i .* t.^2 / 2 + 1i .* 6 .* t) .* exp( -1 .* (t - 4).^2 / 10);

df = 0.01;
f=-5:df:5;

% ============
% Show result
% ============

y = WDF(s + r, t, f);

y = y.';
figure;
image(t, f, abs(y)/max(max(abs(y)))*400);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',12);
ylabel('Frequency (Hz)','Fontsize',12);
xlabel('Time (Sec)','Fontsize',12);
title('STFT of x(t)','Fontsize',12);


% ===
% WDF
% ===
function y = WDF(x, t, f)

    T = length(t);  F = length(f);  % T,F
    dt = t(2) - t(1);   % dt
    df = f(2) - f(1);   % df
    
    mo = round(f/df);  % n_o, m_o

    N = round(1/dt/df/2);
    y = zeros(T, F); 

    f = find(x);
    n1 = f(1);
    n2 = f(end);

    for n = 1:T
        Q = min(n2 - n, n - n1);
        c1 = zeros(1, N);

        for q = 0: 2 * Q
            c1(q + 1) = x(n + q - Q) * conj(x(n - q + Q));
        end

        C = fft(c1);

        for m = 1:F
            y(n, m) = C(mod(mo(m), N)+1)*exp(-1i * 2 * pi * Q * mo(m)/N)* 2 *dt;
        end
    end

end


