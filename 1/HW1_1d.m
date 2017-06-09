clear all; close all; clc

fs = 3*10^4;
T = 0.01;
N = fs*T;
over_smpl = 2;
% M = over_smpl*N;
M = 2^nextpow2(over_smpl*N); % FFT size (zero-padding & over_smpl)
t = (1:N)'/fs;

v = 200/3.6;
L = 200;
fo = 10^4;
Fo = 10^10;
c = 3*10^8;
df = 2*v/c*Fo;
snr = (-20:3:30); % SNR in dB scale
tol = 5/3.6; % Accepted tolerance in km/h
CRB = 12/(N*(N^2 - 1)).*(10.^(-snr/10))*(0.5*c/Fo*fs)^2/(4*pi^2);

Nrun = 100;% # simulation run
f_est1 = zeros(Nrun, length(v));
f_est2 = f_est1;
v1 = f_est1;
v2 = f_est1;
MSE1 = zeros(length(snr), length(v));
MSE2 = MSE1;
MSE1v = zeros(length(snr), length(v));
MSE2v = MSE1;
sigma_w = 10^-10;

for iv = 1:length(v)
    for isnr = 1:length(snr)
        A = sqrt(2*10^(snr(isnr)/10))*sqrt(sigma_w);
        for run = 1:Nrun
            % Signal generation
            phi = 2*pi*rand;
            w = 10^-5*randn(N, 1);
            x = A*cos(2*pi*(fo + df(iv))*t + phi) + w;
            
            % Freq. estimation: DTF on M samples
            S = (abs(fft(x,M)).^2)/M; % DFT on M>N samples (zero-padding)
            [peak, peak_pos] = findpeaks(S(2:end/2), 'SORTSTR', 'descend'); %,'NPeaks',2,'MinPeakDistance',width(iv));
            f_est1(run, iv) = peak_pos(1);
            v1(run, iv)= ((f_est1(run, iv)/M)*fs - fo)*c/(2*Fo);
            
            % Freq. estimation: quadratic interpolation
            f_cent = f_est1(run, iv) + 1;
            Num = S(f_cent - 1) - S(f_cent + 1);
            Den = S(f_cent - 1) + S(f_cent + 1) - 2*S(f_cent);
            f_est2(run, iv) = f_cent + 0.5*Num/Den - 1;
            v2(run, iv) = ((f_est2(run, iv)/M)*fs - fo)*c/(2*Fo);
        end
        MSE1v(isnr, iv) = mean((v1(:, iv) - v(iv)).^2);
        MSE2v(isnr, iv) = mean((v2(:, iv) - v(iv)).^2);
        %         if (MSE2v(isnr, iv) >= tol && MSE2v(isnr, iv) >= 0.05*v(iv))
        %             MSE2v(isnr, iv) = 100000000;
        %         end
    end
end

upper = 10*log10(tol)*ones(1, length(snr));
lower = 0.05*v*ones(1, length(snr));

%figure(1)
for iv = 1:length(v)
    figure()
    semilogy(snr, MSE1v, '-', snr, MSE2v, '-*', snr, CRB, '--', ...
        snr, lower, '-', snr, lower(iv), '-k')
    xlabel('SNR [dB]')
    ylabel('MSE1_v')
    title([num2str(3.6*v(iv)) ' km/h speed'])
    legend('MSE', 'MSE with quad. interpol', 'Cramer Rao Bound', 'Upper limit', 'Lower limit')
end