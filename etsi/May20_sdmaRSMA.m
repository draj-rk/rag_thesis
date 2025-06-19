%%
% Channels

h1 = [1 1 1 1]';

theta1 = 1*pi/9; % [1, 2, 3]
gamma = 1;

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
% Sum-rate vs. SNR

SNR = [0:5:20];
sdma_mrt = zeros(1, length(SNR));
rs_mrt = zeros(1, length(SNR));

for s = 1:length(SNR)
    
    snr_dB = SNR(s);
    snr_ = 10^(snr_dB/10);

    sdma_mrt(s) = sum_rate_sdma(H, w_mrt, snr_);
    rs_mrt(s) = sum_rate_rsma(H, w_mrt, snr_);
end

plot(SNR, sdma_mrt,'-o');
hold on
plot(SNR, rs_mrt,'-o');
xlabel('SNR (dB)')
ylabel('Sum Rate (bits/s/Hz)')
legend('SDMA - MRT', 'RSMA - MRT')
grid on