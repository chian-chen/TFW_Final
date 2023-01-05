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

y = Ambiguity(s + r, t, f);

y = y.';
figure;
image(t, f, abs(y)/max(max(abs(y)))*400);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',12);
ylabel('Frequency (Hz)','Fontsize',12);
xlabel('Time (Sec)','Fontsize',12);
title('STFT of x(t)','Fontsize',12);

% ===============
% low pass filter


[m, n] = size(y);
Filter = zeros(m, n);
centerf = m / 2;
centert = n / 2;
alpha = 1;

for i = 1:m
    for j = 1:n
        Filter(i, j) = sinc((i - centerf)*df * (j - centert)*dt)*exp(-2 * pi * alpha * ((j - centert)*dt)^2); 
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

% =======
% To WDF
% =======

% =================================================
% Find the B and C in p.185 B for eta and C for tau
% =================================================

% =======
% Failed
% =======

