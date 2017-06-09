clear all;% close all; clc
format long

% Sampling param
fs = 4e4;
T = 0.0165;
N = fs*T;
over_smpl = 2;
M = 2^nextpow2(over_smpl*N); % FFT size (zero-padding & over_smpl)
t = (1:N)'/fs;

% Param
v = [10 50 100 200]/3.6;
fo = 1e4;
Fo = 1e10;
c = 3e8;
df = 2*v/c*Fo;
sigma_w = 10^-10;
width = df*0.5*M/fs;
snr = -20:3:30; % SNR in dB scale
CRB = 12/(N*(N^2 - 1)).*(10.^(-snr/10))*(0.5*c/Fo*fs)^2/(4*pi^2);

Nrun = 100;

% Preallocate
f_est1 = zeros(1,Nrun);
f_est2 = f_est1;
v1 = f_est1;
v2 = f_est1;
MSE1 = zeros(length(snr), length(v));
MSE2 = MSE1;
MSE1v = MSE1;
MSE2v = MSE1;

for iv = 1:length(v)
    for isnr = 1:length(snr)
        A1 = sqrt(2*10^(snr(isnr)/10))*sqrt(sigma_w);
        A2 = 1*A1;
        for run = 1:Nrun
            % Signal generation
            phi1 = 2*pi*rand;
            phi2 = 2*pi*rand;
            w = 10^-5*randn(N,1);
            x = A1*cos(2*pi*(fo + df(iv))*t + phi1) + A2*cos(2*pi*fo*t + phi2) + w;
            x = x.*blackman(N);
            
            % Freq. estimation: DTF on M samples
            S = (abs(fft(x, M)).^2)/N; % DFT on M>N samples (zero-padding)
            [peak, peak_pos] = findpeaks(S(2:end/2), 'SORTSTR', 'descend'); %, 'MinPeakDistance', width(iv));
            if peak_pos(1) < peak_pos(2)
                f_est1(run) = (peak_pos(2));
            else
                f_est1(run) = (peak_pos(1));
            end
            v1(run) = ((f_est1(run)/M)*fs - fo)*c/(2*Fo);
            
            % Freq. estimation: quadratic interpolation
            f_cent = f_est1(run) + 1;
            Num = S(f_cent - 1) - S(f_cent + 1);
            Den = S(f_cent - 1) + S(f_cent + 1) - 2*S(f_cent);
            f_est2(run) = f_cent + 0.5*Num/Den - 1;
            v2(run) = ((f_est2(run)/M)*fs - fo)*c/(2*Fo);
        end
        MSE1(isnr, iv) = mean((f_est1/M - (fo - df(iv))/fs).^2);
        MSE2(isnr, iv) = mean((f_est2/M - (fo - df(iv))/fs).^2);
        MSE1v(isnr, iv) = mean((v1 - v(iv)).^2);
        MSE2v(isnr, iv) = mean((v2 - v(iv)).^2);
    end
end

% Plot results
for iv = 1:length(v)
    figure()
    semilogy(snr, MSE1v(:, iv), '-', snr, MSE2v(:, iv), '-*', ...
        snr, CRB, '--')
    xlabel('SNR1 [dB]')
    ylabel('MSE for velocity')
    title(['MSE estimation vs SNR for ' num2str(3.6*v(iv)) ' km/h speed'])
    legend('MSE', 'MSE with quad. interpol', 'Cramer Rao Bound')
end
