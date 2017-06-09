clear all; %close all; clc

% Sampling param
fs = 4e4;
T = 0.01;
N = fs*T;
over_smpl = 2;
M = 2^nextpow2(over_smpl*N); % FFT size (zero-padding & over_smpl)
t = (1:N)'/fs;

% Param
At = 1;
v = [10 50 100 150]/3.6;
L = 20:5:200;
LDec = 10*log10(L);
fo = 1e4;
Fo = 1e10;
c = 3e8;
df = 2*v/c*Fo;
refl = 1/20;
sigma_w = 10^-10;

Nrun = 1000;

% Preallocate
f_est1 = zeros(1, Nrun);
f_est2 = f_est1;
v1 = f_est1;
v2 = f_est1;
snr = zeros(1, length(L));
MSE1 = zeros(length(L), length(v));
MSE2 = MSE1;
MSE1v = MSE1;
MSE2v = MSE1;
A1 = MSE1;

for iv = 1:length(v)
    for il = 1:length(L)
        snr(il) = 10*log10((At^2*refl^2)/(2*(L(il)^4)*sigma_w));
        A1(il) = sqrt(2*10^(snr(il)/10)*(sigma_w));
        A2 = 1*A1; % 1 or 10
        for run = 1:Nrun
            % Signal generation
            phi1 = 2*pi*rand;
            phi2 = 2*pi*rand;
            w = sqrt(sigma_w)*randn(N, 1);
            x = A1(il)*cos(2*pi*(fo + df(iv))*t + phi1) + A2(il)*cos(2*pi*fo*t + phi2) + w;
%             x = x.*blackman(N);
            
            % Freq. estimation: DTF on M samples
            S = (abs(fft(x, M)).^2)/N; %DFT on M>N samples (zero-padding)
            [peak, peak_pos] = findpeaks(S(2:end/2), 'SORTSTR', 'descend');
            if peak_pos(1) < peak_pos(2)
                f_est1(run) = (peak_pos(2));
            else
                f_est1(run) = (peak_pos(1));
            end
            v1(run) = ((f_est1(run)/M)*fs-fo)*c/(2*Fo);
            
            % Freq. estimation: quadratic interpolation
            f_cent = f_est1(run) + 1;
            Num = S(f_cent - 1)-S(f_cent + 1);
            Den = S(f_cent - 1) + S(f_cent + 1)-2*S(f_cent);
            f_est2(run) = f_cent + 0.5*Num/Den - 1;
            v2(run) = ((f_est2(run)/M)*fs - fo)*c/(2*Fo);
        end
        MSE1(il, iv) = mean((f_est1/M - (fo - df(iv))/fs).^2);
        MSE2(il, iv) = mean((f_est2/M - (fo - df(iv))/fs).^2);
        MSE1v(il, iv) = mean((v1 - v(iv)).^2);
        MSE2v(il, iv) = mean((v2 - v(iv)).^2);
    end
end

CRB = 12/(N*(N^2 - 1)).*(10.^(-snr/10))*(0.5*c/Fo*fs)^2/(4*pi^2);

for iv = 1:length(v)
    figure()
    semilogy(LDec, MSE1v(:, iv), '-', LDec, MSE2v(:, iv), '-*', LDec, CRB, '--');
    xlabel('Length [dB]')
    ylabel('MSE for velocity')
    title(['MSE estimation vs SNR for ' num2str(3.6*v(iv)) ' km/h speed'])
end