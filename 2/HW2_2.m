%% 2a)
clear all; close all; clc;
x = importdata('../../Hw2.mat');

% Param
N = length(x);
M = 1;
Ntot = N*M;
K = 20;

% Windowing
win = bartlett(N);
win = sqrt(N)*win./sqrt(sum(win.^2));

% Preallocate
X = zeros(N, M);

for k = 1:K
    % Calculate periodogram
    for m = 1:M
        X(:, m) = fft(win.*x(1+N*(m-1):N+N*(m-1), k), N)';
        X(:, m) = fftshift(X(:, m));
    end
    S_est = (1/N)*(abs(X).^2);
    
    % Estimate frequencies
    [peak, peak_pos] = findpeaks((S_est), 'SORTSTR', 'descend');
    Pw = sum(abs(x(:, k)).^2);    % Total power
    idx = find(10*log10(peak/Pw) > -20);      % Indeces of relevant freqs
    peak_pos = sort(peak_pos(idx));
    freq(k, :) = (peak_pos - 1 - Ntot*0.5)./Ntot;     % Normalized frequencies
    amp(k, :) = sqrt(S_est(peak_pos));
end

% Plot frequencies
figure()
plot(freq, 'o')
xlabel('Iteration')
ylabel('Normalized Frequency')
legend('signal1', 'signal2', 'signal3')
% From the last figure we recognize the same 3 sinusoids for all iteration with
% different amplitudes due to noise realization. Thus, the frequencies are
% any row vector of Target_freq (Normalized)

% Plot amplitude
figure()
plot(1:K, amp)
xlabel('Iteration')
ylabel('Amplitude')
legend('signal1', 'signal2', 'signal3')

%% 2b)

% Param
N = length(x);
M = 1;
Ntot = N*M;
K = 80;

for k = 21:K
    % Calculate periodogram
    for m = 1:M
        X(:, m) = fft(win.*x(1+N*(m-1):N+N*(m-1), k), Ntot);
        X(:, m) = fftshift(X(:, m));
    end
    S_est = (1/N)*(abs(X).^2);
    
    % Estimate frequencies
    [~, peak_pos] = findpeaks((S_est), 'SORTSTR', 'descend');
    peak_pos = sort(peak_pos(idx));
    freq(k, :) = (peak_pos - 1 - Ntot*0.5)./Ntot;       % Normalized frequencies
    amp(k, :) = sqrt(S_est(peak_pos));
end

% for k = 21:K
%     stddev1 = zeros(3, 1);
%     stddev2 = zeros(2, 1);
%     i_freq = 1:3
%     for i = i_freq
%         stddev1(i) = std([freq(21:(k-1), 1); freq(k, i)])
%     end
%     [~, freq_idx] = min(stddev1)
%     freq(k, 1) = freq(k, i_freq(freq_idx));
% 
%     i_freq(freq_idx) = []
%     for i = 1:length(i_freq)
%         stddev2(i) = std([freq(21:(k-1), 2); freq(k, i_freq(i))])
%     end
%     [~, freq_idx] = min(stddev2)
%     freq(k, 2) = freq(k, i_freq(freq_idx));
% 
%     i_freq(freq_idx) = []
%     freq(k, 3) = freq(k, freq_idx);
% end

% Switch two signals to better match behaviour
freq_tmp = freq(51:K, 1);
freq(51:K, 1) = freq(51:K, 2);
freq(51:K, 2) = freq_tmp;

% Plot freq
figure()
plot(21:K, freq(21:K, :))
axis([21 80 -0.5 0.5])
title('Sinusoids detected Vs. iteration')
xlabel('Iteration')
ylabel('Normalized Frequency')
legend('signal1', 'signal2', 'signal3')
