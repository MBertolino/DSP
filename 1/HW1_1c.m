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

% 1c)
v = [1 10 50 100 150]/3.6;
L = 20:5:200;
L_dec = 10*log10(L);
At = 1;
refl = 1/20;

% Param
fo = 1e4; % Tx freq
Fo = 1e10; % Dopplerfreq
c = 3e8;
df = 2*v/c*Fo; % Doppler contribution
sigma_w = 10^-10;

Nrun = 1000; % MC-runs

% Preallocate
Ar = zeros(length(L), 1);
f_est1 = zeros(1, Nrun);
f_est2 = f_est1;
v1 = f_est1;
v2 = f_est1;
MSE1 = zeros(length(L), length(v));
MSE2 = MSE1;
MSE1v = MSE1;
MSE2v = MSE1;
snr = zeros(length(L), 1);

for iv = 1:length(v)
    for il = 1:length(L)
        snr(il) = 10*log10((At^2*refl^2)/(2*(L(il)^4)*sigma_w));
        Ar(il) = sqrt(2*10^(snr(il)/10)*(sigma_w));
        
        for run = 1:Nrun
            % Signal generation
            phi = 2*pi*rand;
            w = 10^-5*randn(N, 1);
            x = Ar(il)*cos(2*pi*(fo + df(iv))*t + phi) + w;
            
            % Freq. estimation: DTF on M samples
            S = (abs(fft(x, M)).^2)/M; % DFT on M>N samples (zero-padding)
            [~, peak_pos] = findpeaks(S(2:end/2), 'SORTSTR', 'descend');
            f_est1(run) = peak_pos(1);
            v1(run) = ((f_est1(run)/M)*fs - fo)*c/(2*Fo);
            
            % Use quad interpolation to improve estimate
            f_cent = f_est1(run) + 1;
            Num = S(f_cent-1) - S(f_cent+1);
            Den = S(f_cent-1) + S(f_cent+1) - 2*S(f_cent);
            f_est2(run) = f_cent + 0.5*Num/Den - 1;
            v2(run) = ((f_est2(run)/M)*fs - fo)*c/(2*Fo);
        end
        
        % Calculate MSE
        MSE1(il, iv) = mean((f_est1/M - (fo - df(iv))/fs).^2);
        MSE2(il, iv) = mean((f_est2/M - (fo - df(iv))/fs).^2);
        MSE1v(il, iv) = mean((v1 - v(iv)).^2);
        MSE2v(il, iv) = mean((v2 - v(iv)).^2);
    end
end

CRB = 12/(N*(N^2 - 1)).*(10.^(-snr/10))*(0.5*c/Fo*fs)^2/(4*pi^2);

% Plot results
% 1c) MSE vs length L
for iv = 1:length(v)
    figure()
    semilogy(L_dec, MSE1v(:, iv), '-', L_dec, MSE2v(:, iv), '-*', ...
        L_dec, CRB, '--')
    xlabel('Distance [log(m)]')
    ylabel('MSE for velocity')
    title(['MSE vs distance for ' num2str(3.6*v(iv)) ' km/h speed'])
    legend('MSE', 'MSE with quad. interpol', 'Cramer Rao Bound')
end
