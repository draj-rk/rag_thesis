%%
% Channels

h1 = 1+0*1i * [1 1 1 1]';

theta1 = 1*pi/9; % [1, 2, 3]
gamma = 1;

thetaN = 1;

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
% ZF preoders

w_zf = H_norm*inv(H_norm'*H_norm);
w_zf = w_zf ./ vecnorm(w_zf); % Normalize wrt each column
%%
% Sum-rate vs. SNR

SNR = [0:1:20];

sdma_mrt = zeros(1, length(SNR));
rs_mrt = zeros(1, length(SNR));

sdma_zf = zeros(1, length(SNR));
rs_zf = zeros(1, length(SNR));

for s = 1:length(SNR)
    
    snr_dB = SNR(s);
    snr_ = 10^(snr_dB/10);
    
    % MRT precoders
    sdma_mrt(s) = sum_rate_sdma(H, w_mrt, snr_);
    rs_mrt(s) = sum_rate_rsma(H, w_mrt, snr_);
    
    % ZF precoders
    sdma_zf(s) = sum_rate_sdma(H, w_zf, snr_);
    rs_zf(s) = sum_rate_rsma(H, w_zf, snr_);
end

figure
plot(SNR, sdma_mrt,'ro-');
hold on
plot(SNR, rs_mrt,'bo-');
hold on
plot(SNR, sdma_zf,'r^-');
hold on
plot(SNR, rs_zf,'b^-');

xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
legend('SDMA - MRT', 'RSMA - MRT', 'SDMA - ZF', 'RSMA - ZF')
grid on
title('Sum Rate vs. SNR')

%% 
% Gain of RSMA over SDMA

gain_mrt = (rs_mrt - 1*sdma_mrt)./sdma_mrt*100;
gain_zf = (rs_zf - 1*sdma_zf)./sdma_zf*100;
%%
figure
plot(SNR, gain_mrt,'ro-');
hold on
plot(SNR, gain_zf,'bo-');
xlabel('SNR (dB)')
ylabel('Sum Rate Gain')
legend('MRT', 'ZF')
grid on
title(sprintf('Sum Rate Gain \\gamma = %d, \\theta = %d\\pi/9', gamma, thetaN))
