%%
clc
clear all

%%
% Channels
h1 = 1+0*1i * [1 1 1 1]';

gamma = 0.3;
thetaN = 3;

theta1 = thetaN*pi/9; % [1, 3]

h2 = gamma*[1 exp(1i*theta1) exp(1i*2*theta1) exp(1i*3*theta1)]';

H = [h1, h2];
%%
% Normalized H

h1_norm = h1 / norm(h1);
h2_norm = h2 / norm(h2);

H_norm = [h1_norm h2_norm];

%%
% MRT precoders

w_mrt = H_norm;

%%
% ZF precoders

w_zf = H_norm*inv(H_norm'*H_norm);
w_zf = w_zf ./ vecnorm(w_zf); % Normalize wrt each column
%%
% Sum-rate vs. SNR

SNR = [0:1:25];

sdma_mrt = zeros(1, length(SNR));
sdma_mrt_max = zeros(1, length(SNR));
sdma_mrt_min = zeros(1, length(SNR));
rs_mrt = zeros(1, length(SNR));
rs_mrt_max = zeros(1, length(SNR));
rs_mrt_min = zeros(1, length(SNR));
oma_mrt = zeros(1, length(SNR));
oma_mrt_max = zeros(1, length(SNR));
oma_mrt_min = zeros(1, length(SNR));
noma_mrt = zeros(1, length(SNR));
noma_mrt_max = zeros(1, length(SNR));
noma_mrt_min = zeros(1, length(SNR));

sdma_zf = zeros(1, length(SNR));
sdma_zf_max = zeros(1, length(SNR));
sdma_zf_min = zeros(1, length(SNR));
rs_zf = zeros(1, length(SNR));
rs_zf_max = zeros(1, length(SNR));
rs_zf_min = zeros(1, length(SNR));

sdma_mm = zeros(1, length(SNR));
sdma_mm_max = zeros(1, length(SNR));
sdma_mm_min = zeros(1, length(SNR));
rs_mm = zeros(1, length(SNR));
rs_mm_max = zeros(1, length(SNR));
rs_mm_min = zeros(1, length(SNR));

for s = 1:length(SNR)
    
    snr_dB = SNR(s);
    snr_ = 10^(snr_dB/10);
    
    % % MRT precoders
    % [sdma_mrt(s), sdma_mrt_max(s), sdma_mrt_min(s)] = sum_rate_sdma(H, w_mrt, snr_);
    % [rs_mrt(s), rs_mrt_max(s), rs_mrt_min(s)] = sum_rate_rsma(H, w_mrt, snr_);
    [oma_mrt(s), oma_mrt_max(s), oma_mrt_min(s)] = sum_rate_oma(H, w_mrt, snr_);
    [noma_mrt(s), noma_mrt_max(s), noma_mrt_min(s)] = sum_rate_NOMA(H, w_mrt, snr_);
    % 
    % % ZF precoders
    % [sdma_zf(s), sdma_zf_max(s), sdma_zf_min(s)] = sum_rate_sdma(H, w_zf, snr_);
    % [rs_zf(s), rs_zf_max(s), rs_zf_min(s)] = sum_rate_rsma(H, w_zf, snr_);

    % MMSE precoders
    reg_ = eye(2)*(1/(snr_/2));
    w_mm = H*inv(H'*H + reg_);
    w_mm = w_mm ./ vecnorm(w_mm); % Normalize wrt each column

    [sdma_mm(s), sdma_mm_max(s), sdma_mm_min(s)] = sum_rate_sdma(H, w_mm, snr_);
    [rs_mm(s), rs_mm_max(s), rs_mm_min(s)] = sum_rate_rsma_mmse(H, H_norm, w_mm, snr_);

end

%%
% Average user rates
figure

plot(SNR, sdma_mm,'r+-');
hold on
plot(SNR, oma_mrt,'mx--');
hold on
plot(SNR, rs_mm,'bo-');
hold on
plot(SNR, noma_mrt,'g^--');

ylim([0 8])
xlabel('SNR (dB)')
ylabel('Average user-rates')
legend('SDMA - MMSE', 'OMA', 'RSMA - MMSE', 'NOMA', 'Location', 'northwest')
grid on
title(['MMF maximisation, $\gamma = ' num2str(gamma) ',\ \theta = ' num2str(thetaN) '\frac{\pi}{9}$'], 'Interpreter', 'latex');
%%
% Weak user rates
figure

plot(SNR, sdma_mm_min,'r+-');
hold on
plot(SNR, oma_mrt_min,'mx--');
hold on
plot(SNR, rs_mm_min,'bo-');
hold on
plot(SNR, noma_mrt_min,'g^--');

ylim([0 8])
xlabel('SNR (dB)')
ylabel('Weak user-rates')
legend('SDMA - MMSE', 'OMA', 'RSMA - MMSE', 'NOMA', 'Location', 'northwest')
grid on
title(['MMF maximisation, $\gamma = ' num2str(gamma) ',\ \theta = ' num2str(thetaN) '\frac{\pi}{9}$'], 'Interpreter', 'latex');
