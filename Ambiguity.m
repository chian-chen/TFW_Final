[a1, fs] = audioread('Chord.wav');
x = a1(:,1).';    % only extract the first channel
tau = 0:1/fs:length(x)/fs;
dt = 0.01; df= 1; sgm= 200;
t= 0:dt:max(tau); f= 20:df:1000;

tic
y = Gabor(x, tau, t, f, sgm);
toc

y = y.';
figure;
image(t, f, abs(y)/max(max(abs(y)))*400);
colormap(gray(256));
set(gca,'Ydir','normal');
set(gca,'Fontsize',12);
ylabel('Frequency (Hz)','Fontsize',12);
xlabel('Time (Sec)','Fontsize',12);
title('STFT of x(t)','Fontsize',12);


function y = Gabor(x, tau, t, f, sgm)
    T = length(t);  F = length(f);  % T,F
    dt = t(2) - t(1);               % dt
    df = f(2) - f(1);               % df
    dtau = tau(2) - tau(1);         % dtau
    S = ceil(dt/dtau);              % S
    
    n_o = t(1)/dt;  m_o = f(1)/df;  % n_o, m_o
    N = 1/dtau/df;  Q = ceil(1.9143/dtau/sqrt(sgm));    % N, Q
    
    y = zeros(T, F); 
    
    % O(T)
    for n = n_o: n_o + T - 1
        % Step: determine x1(q)     O(2Q)
        x1 = zeros(1, N);
        
        for q = 0: 2*Q
            index = n * S + q - Q;
            if index >= 0 && index < length(x)
                x1(q + 1) = x(index + 1) * exp(-sgm*pi*((q-Q)*dtau)^2) * sgm^(1/4);
            end
        end
  
        % Step: FFT     O(NlogN)
        X = fft(x1);
        
        % Step: Convert X(m) to y(T, F)     O(F)
        for m = m_o : m_o + F - 1
            y(n - n_o + 1, m - m_o + 1) = X(mod(m,N)+1)*exp(-1i*2*pi*(Q - n*S)*m)*dtau;
        end
    end
end