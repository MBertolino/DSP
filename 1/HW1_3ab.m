% Homework 3a
% Digital Signal Processing
%
% Mattias Bertolino
clear all; close all; clc

% Sampling param
c = 3e8;
L = 200;
Tg = 1e-8;
Ng = 21;
ts = Tg/Ng;
tauL = 2*L/c;
T = Tg + tauL; % Time window large enough
N = ceil(T/ts);

% Policeman param - first pulse
L1 = 200; % Where the policeman stands
tau1 = 2*L1/c;
N1 = ceil(tau1/ts);

% Pulse init
sigma_w = 1e-10;
snr = -20:3:30;
A = sqrt(3*10.^(snr./10)*(sigma_w));
g = triang(Ng)';

% Preallocate
Nrun = 100;
x1 = zeros(1, N1);
tod1 = zeros(1, Nrun);
MSE_tod = zeros(1, length(snr));

for isnr = 1:length(snr)
    for run = 1:Nrun
        % Pulse generation
        w1 = sqrt(sigma_w)*randn(1, N);
        x1 = A(isnr)*[zeros(1, N1) g zeros(1, N - Ng - N1)] + w1;
        
        % Time of Delay estimation
        [xcor, lag] = xcorr(x1, g);
        [~, idx] = max(xcor.^2);
        tod1(run) = lag(idx)*ts;
    end
    MSE_tod(isnr) = mean((tod1 - tau1).^2);
end

% Cramer-Rao Bound for Time of Delay
CRB_tod = 1./(2400*10.^(snr/10)/Tg^2);

% 3a) Plot MSE for Time of Delay estimation
figure()
semilogy(snr, MSE_tod, snr, CRB_tod, '--');
xlabel('Signal to noise-ratio');
ylabel('MSE of Time of Delay');
title(['MSE vs SNR for ' num2str(L) ' meters'])
legend('MSE', 'Cramer Rao Bound')


%% 3b) MSE for Speed estimation
dt = 0.5;
v = [10 50 100 150];

% Preallocate
x2 = zeros(1, N1);
tod2 = zeros(1, Nrun);
v_est = zeros(1, Nrun);
MSE_speed = zeros(length(v), length(snr));

for iv = 1:length(v)
    % Second pulse
    tau2 = (2*(L1 - dt*v(iv)/3.6)/c);
    N2 = ceil(tau2/ts);
    for isnr = 1:length(snr)
        for run = 1:Nrun
            % Pulse generation
            w1 = sqrt(sigma_w)*randn(N, 1);
            w2 = sqrt(sigma_w)*randn(N, 1);
            x1 = A(isnr)*[zeros(1, N1) g zeros(1, N - Ng - N1)]' + w1;
            x2 = A(isnr)*[zeros(1, N2) g zeros(1, N - Ng - N2)]' + w2;
            
            % Time of Delay estimation 1
            [xcor, lag] = xcorr(x1, g);
            [~, idx] = max(xcor.^2);
            tod1(run) = lag(idx)*ts;
            
            % Time of Delay estimation 2
            [xcor, lag] = xcorr(x2, g);
            [~, idx] = max(xcor.^2);
            tod2(run) = lag(idx)*ts;
            
            % Speed estimation
            v_est(run) = c*(tod1(run) - tod2(run))/dt*0.5*3.6;
        end
        MSE_speed(iv, isnr) = mean((v_est - v(iv)).^2);
    end
end

% CRB for Speed estimation
CRB_speed = 2./(24*10.^(snr/10)/Tg^2)*(c/(2*dt))^2;

% 3b) Plot MSE for Speed estimation
for iv = 1:length(v)
    figure()
    semilogy(snr, MSE_speed(iv, :), snr, CRB_speed);
    xlabel('Signal to noise-ratio');
    ylabel('MSE of Time of Delay');
    title(['MSE vs SNR for ' num2str(v(iv)) ' km/h'])
    legend('MSE', 'Cramer Rao Bound')
end