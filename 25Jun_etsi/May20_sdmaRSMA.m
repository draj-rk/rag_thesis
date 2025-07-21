%%
clc
clear all

%%
% Channels
h1 = 1+0*1i * [1 1 1 1]';

gamma = 0.3;
thetaN = 3;

theta1 = thetaN*pi/9; % [1, 2, 3]

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

SNR = [0:1:20];

sdma_mrt = zeros(1, length(SNR));
sdma_mrt_max = zeros(1, length(SNR));
sdma_mrt_min = zeros(1, length(SNR));
rs_mrt = zeros(1, length(SNR));
rs_mrt_max = zeros(1, length(SNR));
rs_mrt_min = zeros(1, length(SNR));
oma_mrt = zeros(1, length(SNR));
oma_mrt_max = zeros(1, length(SNR));
oma_mrt_min = zeros(1, length(SNR));

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
    
    % MRT precoders
    [sdma_mrt(s), sdma_mrt_max(s), sdma_mrt_min(s)] = sum_rate_sdma(H, w_mrt, snr_);
    [rs_mrt(s), rs_mrt_max(s), rs_mrt_min(s)] = sum_rate_rsma(H, w_mrt, snr_);
    [oma_mrt(s), oma_mrt_max(s), oma_mrt_min(s)] = sum_rate_oma(H, w_mrt, snr_);
    
    % ZF precoders
    [sdma_zf(s), sdma_zf_max(s), sdma_zf_min(s)] = sum_rate_sdma(H, w_zf, snr_);
    [rs_zf(s), rs_zf_max(s), rs_zf_min(s)] = sum_rate_rsma(H, w_zf, snr_);

    % MMSE precoders
    reg_ = eye(2)*(1/(snr_/2));
    w_mm = H*inv(H'*H + reg_);
    w_mm = w_mm ./ vecnorm(w_mm); % Normalize wrt each column

    [sdma_mm(s), sdma_mm_max(s), sdma_mm_min(s)] = sum_rate_sdma(H, w_mm, snr_);
    [rs_mm(s), rs_mm_max(s), rs_mm_min(s)] = sum_rate_rsma_mmse(H, w_mm, snr_);

end

%%
% figure
% plot(SNR, sdma_mrt,'ro-');
% hold on
% plot(SNR, rs_mrt,'bo-');
% hold on
% plot(SNR, sdma_zf,'r^-');
% hold on
% plot(SNR, rs_zf,'b^-');
% hold on
% plot(SNR, sdma_mm,'r--');
% hold on
% plot(SNR, rs_mm,'b--');
% 
% xlabel('SNR (dB)')
% ylabel('Sum Rate (bits/s/Hz)')
% legend('SDMA - MRT', 'RSMA - MRT', 'SDMA - ZF', 'RSMA - ZF', 'SDMA - MMSE', 'RSMA - MMSE')
% grid on
% title('Sum Rate vs. SNR')

%% 
% Gain of RSMA over SDMA

gain_mrt = (rs_mrt - 1*sdma_mrt)./sdma_mrt*100;
gain_mrt_max = (rs_mrt_max - 1*sdma_mrt_max)./sdma_mrt_max*100;
gain_mrt_min = (rs_mrt_min - 1*sdma_mrt_min)./sdma_mrt_min*100;

gain_zf = (rs_zf - 1*sdma_zf)./sdma_zf*100;
gain_zf_max = (rs_zf_max - 1*sdma_zf_max)./sdma_zf_max*100;
gain_zf_min = (rs_zf_min - 1*sdma_zf_min)./sdma_zf_min*100;

gain_mm = (rs_mm - 1*sdma_mm)./sdma_mm*100;
gain_mm_max = (rs_mm_max - 1*sdma_mm_max)./sdma_mm_max*100;
gain_mm_min = (rs_mm_min - 1*sdma_mm_min)./sdma_mm_min*100;
%%
figure
subplot(1,3,1);
plot(SNR, gain_mrt,'ro-');
hold on
plot(SNR, gain_zf,'bo-');
hold on
plot(SNR, gain_mm,'m--');

xlabel('SNR (dB)')
ylabel('Sum Rate Gain')
legend('MRT', 'ZF', 'MMSE')
grid on
title(sprintf('Gain RSMA over SDMA \\gamma = %.1f, \\theta = %d\\pi/9', gamma, thetaN))

subplot(1,3,2);
plot(SNR, gain_mrt_min,'ro-');
hold on
plot(SNR, gain_zf_min,'bo-');
hold on
plot(SNR, gain_mm_min,'m--');

xlabel('SNR (dB)')
ylabel('Weak user-rate gain')
legend('MRT', 'ZF', 'MMSE')
grid on

subplot(1,3,3);
plot(SNR, gain_mrt_max,'ro-');
hold on
plot(SNR, gain_zf_max,'bo-');
hold on
plot(SNR, gain_mm_max,'m--');

xlabel('SNR (dB)')
ylabel('Strong user-rate gain')
legend('MRT', 'ZF', 'MMSE')
grid on

%%
% Gain of SDMA over OMA
gain_mrt = (sdma_mrt - 1*oma_mrt)./oma_mrt*100;
gain_mrt_max = (sdma_mrt_max - 1*oma_mrt_max)./oma_mrt_max*100;
gain_mrt_min = (sdma_mrt_min - 1*oma_mrt_min)./oma_mrt_min*100;

gain_zf = (sdma_zf - 1*oma_mrt)./oma_mrt*100;
gain_zf_max = (sdma_zf_max - 1*oma_mrt_max)./oma_mrt_max*100;
gain_zf_min = (sdma_zf_min - 1*oma_mrt_min)./oma_mrt_min*100;

gain_mm = (sdma_mm - 1*oma_mrt)./oma_mrt*100;
gain_mm_max = (sdma_mm_max - 1*oma_mrt_max)./oma_mrt_max*100;
gain_mm_min = (sdma_mm_min - 1*oma_mrt_min)./oma_mrt_min*100;

%%
figure
subplot(1,3,1);
plot(SNR, gain_mrt,'ro-');
hold on
plot(SNR, gain_zf,'bo-');
hold on
plot(SNR, gain_mm,'m--');

xlabel('SNR (dB)')
ylabel('Sum Rate Gain')
legend('MRT', 'ZF', 'MMSE')
grid on
title(sprintf('Gain of SDMA over OMA \\gamma = %.1f, \\theta = %d\\pi/9', gamma, thetaN))

subplot(1,3,2);
plot(SNR, gain_mrt_min,'ro-');
hold on
plot(SNR, gain_zf_min,'bo-');
hold on
plot(SNR, gain_mm_min,'m--');

xlabel('SNR (dB)')
ylabel('Weak user-rate gain')
legend('MRT', 'ZF', 'MMSE')
grid on

subplot(1,3,3);
plot(SNR, gain_mrt_max,'ro-');
hold on
plot(SNR, gain_zf_max,'bo-');
hold on
plot(SNR, gain_mm_max,'m--');

xlabel('SNR (dB)')
ylabel('Strong user-rate gain')
legend('MRT', 'ZF', 'MMSE')
grid on
