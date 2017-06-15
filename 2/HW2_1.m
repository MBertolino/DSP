%% 1a) Periodogram

clear all; close all; clc;
N = 1e4;
M = 30;
Ntot = M*N;
freq = linspace(0, 1, N/2);

% Window
win = bartlett(N);
win = sqrt(N)*win./sqrt(sum(win.^2));

% Given LP random process
[b, a] = cheby1(15, 0.1, 0.2); % 15th order Chebychev filter
A = 1e-4;
omega_o = pi*0.35;
z = filter(b, a, randn(Ntot, 1));
x = z + A*cos(omega_o*(1:Ntot)');

% Calculate Periodogram
X = zeros(M, N);
for m = 1:M
    X(m, :) = fft(win.*x(1+N*(m-1):N+N*(m-1)), N)';
    X(m, :) = fftshift(X(m, :));
end
S_est = 10*log10((1/N)*mean(abs(X).^2));

% Plot results
figure()
plot(freq, S_est(end/2 + 1:end))
xlabel('Normalized frequency (\times\pi rad/sample)')
ylabel('Power/frequency (dB/rad/sample)')
title('Power Spectral Density Estimate')
grid on;
hold on;
% figure()
% periodogram(x)
% figure()
% pwelch(x)

%% 1b) AR(Na) Autoregressive
Na = 10:10:40;

% Method 1 (CHEATING)
% figure()
% pyulear(x, Na(1), N)
% title('AR cheat')

% Preallocate
S = zeros(1, length(freq));

for i = 1:length(Na)
    % Method 2 (w/ Levinson regression)
    % Rxx = xcorr(x, 'coeff');
    % cent = round((length(Rxx)+1)/2);
    % Rxx(1:Na(i)+1) = Rxx(cent:cent+Na(i));
    % p = levinson(Rxx, Na(i));
    % [H, w] = freqz(1, p, N);
    % % figure()
    % plot (w/pi, 10*log10(H));
    % title('AR with levinson');
    
    % Method 3 (just plain)
    % Setup the Yule-Walker equations
    r = xcorr(x, Na(i), 'coeff');
    %     cent = round((length(r)+1)/2);
    %     r(1:Na(i)+1) = r(cent:cent+Na(i));
    r(1:Na(i)+1) = r(Na(i)+1:end);
    R = toeplitz(r(1:Na(i)));
    r = r(2:Na(i)+1);
    p2 = [1; -R\r];
        
    % Estimate AR
    for j = 1:length(freq)
        z = exp(-1i*pi*freq(j));
        k = 0:Na(i);
        den = sum(p2'.*z.^(-k));
        
        S(j) = abs(1/den)^2;
    end
    S = S./(abs(S(1)));
    
    % Alternative way to calculate
    %     coef = @(z) [];
    %     for k = 0:Na(i)
    %         coef = @(z) [coef(z) z.^(-k)];
    %     end
    %     A = @(z) 1/sum(coef(z)*p2);
    %     S = @(z) A(z)*conj(A(1/z));
    %
    %     s1 = zeros(1, length(freq));
    %     for j = 1:length(freq)
    %         s1(j) = S(exp(-1i*pi*freq(j)));
    %     end
    
    % Plot AR
    plot(freq, 10*log10(S))
end
legend('Periodogram', 'AR(10)', 'AR(20)', 'AR(30)', 'AR(40)')

% Plot Pole-Zero Map
iir = tf(1, p2');
figure()
pzmap(iir)

%% 1c) MA(Nb) Moving Average

Nb = 100;
h = boxcar(Nb);
% ma = conv(yup,x);

% plot(ma)

m = length(x); % Length of input signal
n = length(h); % Length of window
X = [fft(x)', zeros(1, Nb)]; % Make x and h same length vector X, H
H = [h', zeros(1, m)];
for i = 1:Nb+Ntot-1 % Create output vector
    if i == 1
        Y(1) = (X(1)*H(1))/Nb;
    elseif i <= Nb
        Y(i) = Y(i-1)+X(i)/Nb;
    else
        Y(i) = Y(i-1)+X(i)/Nb - X(i-Nb)/Nb;
    end
end