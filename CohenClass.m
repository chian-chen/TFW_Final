% ===================================
% Two function s and r is from p.96
% ===================================

dt = 0.01;
df = 0.01;
f=-5:df:5;

t = -9:dt:1;
s = exp(1i * t.^2 / 10 - 1i * 3 * t);
s = [zeros(1, 100) s zeros(1, 900)];

t = -10:dt:10;
r = exp(1i .* t.^2 / 2 + 1i .* 6 .* t) .* exp( -1 .* (t - 4).^2 / 10);


% ===================
% Ambiguity Function
% ===================

y = Ambiguity(s, t, f);

y = y.';
figure;
image(t, f, abs(y)/max(max(abs(y)))*400);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',12);
ylabel('Frequency (Hz)','Fontsize',12);
xlabel('Time (Sec)','Fontsize',12);
title('STFT of x(t)','Fontsize',12);


% low pass filter

[m, n] = size(y);
Filter = zeros(m, n);
centerf = m / 2;
centert = n / 2;

for i = 1:m
    for j = 1:n
        Filter(i, j) = sinc((i - centerf)*df * (j - centert)*dt)*exp(-2 * pi * ((j - centert)*dt)^2); 
    end
end

figure;
image(t, f, abs(Filter)/max(max(abs(Filter)))*400);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',12);
ylabel('Frequency (Hz)','Fontsize',12);
xlabel('Time (Sec)','Fontsize',12);
title('STFT of x(t)','Fontsize',12);

y = y .* Filter;

% ================
% To WDF (Failed)
% ================




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
        
        pp = p;
        if(pp > 1001)
            pp = pp - 1001;
        else
            pp = 1001 - pp;
        end

        Q = n2 - n1 - 2 * pp;
        c1 = zeros(1, N);
        
        for n = 0 : Q
            c1(n + 1) = x(n + n1 + 2 * pp) * conj(x(n + n1));
        end

        C = fft(c1);

        for m = 1:F
            y(p, m) = C(mod(mo(m), N)+1)*exp(-1i * 2 * pi * Q * mo(m)/N)*dt;
        end
    end
end
