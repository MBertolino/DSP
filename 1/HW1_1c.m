% Homework 1a, 1b
% Digital Signal Processing
%
% Mattias Bertolino
clear all; close all; clc

% Sampling param
fs = 4e4; % sampl
T = 0.01;
N = fs*T;
over_smpl = 2;
M = 2^nextpow2(over_smpl*N); % FFT size (zero-padding & over_smpl)
t = (1:N)'/fs;

L = [20 100 200];
% 1a)
v = [0 1 10 50 100 200]/3.6; % km/h to m/s
% 1b)
% v = 100/3.6;

% Param
fo = 1.e4; % Tx freq
Fo = 1e10; % Dopplerfreq
c = 3e8;
df = 2*v/c*Fo; % Doppler contribution
snr = -20:3:30; % SNR in dB scale
sigma_w = 10^-10;

Nrun = 100; % MC-runs

% Preallocate
f_est1 = zeros(1, Nrun);
f_est2 = f_est1;
v1 = f_est1;
v2 = f_est1;
A_power = zeros(length(snr), length(L));
MSE1 = zeros(length(snr), length(v));
MSE2 = MSE1;
MSE1v = MSE1;
MSE2v = MSE1;
CRB = 12/(N*(N^2 - 1)).*(10.^(-snr/10))*(0.5*c/Fo*fs)^2/(4*pi^2);
MSE_floor = (1/3)*(.5/M)^2; % quantization error +/-(.5/M)

for iv = 1:length(v)
    for isnr = 1:length(snr)
        A = sqrt(2*10^(snr(isnr)/10))*sqrt(sigma_w); % Ampl depends on snr
        for run = 1:Nrun
            % Signal generation
            phi = 2*pi*rand;
            w = 10^-5*randn(N, 1);
            x = A*cos(2*pi*(fo + df(iv))*t + phi) + w;
            
            % Freq. estimation: DTF on M samples
            S = (abs(fft(x, M)).^2)/M; % DFT on M>N samples (zero-padding)
            
            [~, peak_pos] = findpeaks(S(2:end/2), 'SORTSTR','descend');
            f_est1(run) = peak_pos(1);
            v1(run) = ((f_est1(run)/M)*fs - fo)*c/(2*Fo);
            
            % Use quad interpolation to improve estimate
            f_cent = f_est1(run) + 1;
            Num = S(f_cent-1) - S(f_cent+1);
            Den = S(f_cent-1) + S(f_cent+1) - 2*S(f_cent);
            f_est2(run) = f_cent + 0.5*Num/Den - 1;
            v2(run) = ((f_est2(run)/M)*fs-fo)*c/(2*Fo);
        end
        
        % Calculate MSE
        MSE1(isnr, iv) = mean((f_est1/M - (fo - df(iv))/fs).^2);
        MSE2(isnr, iv) = mean((f_est2/M - (fo - df(iv))/fs).^2);
        MSE1v(isnr, iv) = mean((v1 - v(iv)).^2);
        MSE2v(isnr, iv) = mean((v2 - v(iv)).^2);
        
        % 1b) MSE vs At^2/2 (Power)
        A_power(isnr, :) = 10*log10(0.5*(20*L.^2*A).^2);
    end
end


% Plot results
% 1a) MSE vs SNR
for iv = 1:length(v)
    % figure()
    % semilogy(snr, MSE1v(:, iv), '-', snr, MSE2v(:, iv), '-*', ...
    %     snr, CRB, '--', snr, MSE_floor*(1+0*snr), ':', ...
    %     snr, ((1/4)/12)*(1+0*snr), ':')
    % xlabel('SNR [dB]')
    % ylabel('MSE for velocity')
    % title(['MSE vs SNR for ' num2str(3.6*v(iv)) ' km/h speed'])
    % legend('MSE', 'MSE with quad. interpol', 'Cramer Rao Bound')
end
