% =========================================
% function from p.251 (to check the result)
% =========================================

t = -5:0.1:5;
x = [zeros(1, 20) ones(1, 61) zeros(1,20)];

figure;
plot(t, x);
ylim([-1, 2])


a = 0.1;
y = frft(x, a);

figure;
plot(t, y);
ylim([-1, 2])


% function y = FRFT(x, a)
% % phi = 0.5 * a * pi
%     a = mod(a,4);
%     y = x;
%     
%     
%     
%     
% end

