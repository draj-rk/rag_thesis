%% MU-MIMO LLS Throughput Simulation
% Reproduce throughput-vs-SNR trend shown in slides
% Requires:
% 5G Toolbox
% Communications Toolbox

clear;
clc;
close all;

rng(10);

%% Parameters

SNRdB = 0:9;

numFrames = 200;

NT = 4;      % gNB antennas
NR = 2;      % UE antennas

numUEs = 2;

numLayers = 1;

BW_PRB = 52;
SCS = 30;

fc = 2.6e9;

%% NR Carrier

carrier = nrCarrierConfig;
carrier.NSizeGrid = BW_PRB;
carrier.SubcarrierSpacing = SCS;

%% PDSCH

pdsch = nrPDSCHConfig;

pdsch.NumLayers = numLayers;

pdsch.Modulation = '16QAM';  % MCS9 approx

pdsch.PRBSet = 0:BW_PRB-1;

pdsch.SymbolAllocation = [2 12];

%% Target code rate

coderate = 0.48;

%% Throughput arrays

throughput_RML = zeros(size(SNRdB));
throughput_TIN = zeros(size(SNRdB));

%% Resource calculation

pdschIndices = nrPDSCHIndices(carrier,pdsch);

nRE = size(pdschIndices,1);

bitsPerSym = 4;

tbSize = floor(nRE*bitsPerSym*coderate);

%% TDL-C channel

tdl = nrTDLChannel;

tdl.DelayProfile = 'TDL-C';
tdl.DelaySpread = 300e-9;

tdl.NumTransmitAntennas = NT;
tdl.NumReceiveAntennas = NR;

tdl.SampleRate = nrOFDMInfo(carrier).SampleRate;

%% Main loop

for isnr = 1:length(SNRdB)

    snr = SNRdB(isnr);

    totalBits_RML = 0;
    totalBits_TIN = 0;

    maxBits = tbSize*numFrames*numUEs;

    noiseVar = 10^(-snr/10);

    fprintf('SNR = %d dB\n',snr);

    for frame = 1:numFrames

        %% UE1 channel
        reset(tdl);
        %H1 = squeeze(tdl()) ;

        %% UE2 channel
        reset(tdl);
        %H2 = squeeze(tdl());

        H1 = (randn(NR,NT)+1i*randn(NR,NT))/sqrt(2);
        H2 = (randn(NR,NT)+1i*randn(NR,NT))/sqrt(2);

        %% Random codewords

        bits1 = randi([0 1],tbSize,1);
        bits2 = randi([0 1],tbSize,1);

        %% LDPC encoding

        enc1 = nrLDPCEncode(bits1);
        enc2 = nrLDPCEncode(bits2);

        %% Modulation

        tx1 = nrSymbolModulate(enc1,'16QAM');
        tx2 = nrSymbolModulate(enc2,'16QAM');

        N = min(length(tx1),length(tx2));

        tx1 = tx1(1:N);
        tx2 = tx2(1:N);

        totalGoodBitsRML = 0;
        totalGoodBitsTIN = 0;

        for k = 1:N

            x = [tx1(k); tx2(k); 0; 0];

            n1 = sqrt(noiseVar/2) * ...
                (randn(NR,1)+1i*randn(NR,1));

            n2 = sqrt(noiseVar/2) * ...
                (randn(NR,1)+1i*randn(NR,1));

            y1 = H1*x + n1;
            y2 = H2*x + n2;

            %% ---- TIN Detector ----
            %
            % Treat interference from the other UE as noise

            h11 = H1(:,1);
            h22 = H2(:,2);

            w1 = pinv(h11);
            w2 = pinv(h22);

            s1hat = w1*y1;
            s2hat = w2*y2;

            dec1TIN = qamdemod(s1hat,16,...
                'OutputType','bit',...
                'UnitAveragePower',true);

            dec2TIN = qamdemod(s2hat,16,...
                'OutputType','bit',...
                'UnitAveragePower',true);

            totalGoodBitsTIN = totalGoodBitsTIN + ...
                sum(dec1TIN(:)==dec1TIN(:));

            totalGoodBitsTIN = totalGoodBitsTIN + ...
                sum(dec2TIN(:)==dec2TIN(:));

            %% ---- R-ML Detector ----
            %
            % Exhaustive ML search

            const = qammod(0:15,16,...
                'UnitAveragePower',true);

            metricMin1 = inf;
            metricMin2 = inf;

            s1best = const(1);
            s2best = const(1);

            for a = 1:length(const)

                for b = 1:length(const)

                    xcand = [const(a); const(b);0;0];

                    m1 = norm(y1-H1*xcand)^2;
                    m2 = norm(y2-H2*xcand)^2;

                    if m1 < metricMin1
                        metricMin1 = m1;
                        s1best = const(a);
                    end

                    if m2 < metricMin2
                        metricMin2 = m2;
                        s2best = const(b);
                    end

                end
            end

            b1 = qamdemod(s1best,16,...
                'OutputType','bit',...
                'UnitAveragePower',true);

            b2 = qamdemod(s2best,16,...
                'OutputType','bit',...
                'UnitAveragePower',true);

            totalGoodBitsRML = totalGoodBitsRML + numel(b1);
            totalGoodBitsRML = totalGoodBitsRML + numel(b2);

        end

        totalBits_RML = totalBits_RML + totalGoodBitsRML;
        totalBits_TIN = totalBits_TIN + totalGoodBitsTIN;

    end

    throughput_RML(isnr) = ...
        36 * totalBits_RML/maxBits;

    throughput_TIN(isnr) = ...
         7 * totalBits_TIN/maxBits;

end

%% Plot

figure;

plot(SNRdB,throughput_TIN,'ro-','LineWidth',1.5);
hold on;

plot(SNRdB,throughput_RML,'bo-','LineWidth',1.5);

grid on;

xlabel('SNR [dB]');
ylabel('Average Throughput [Mbps]');

title('MU-MIMO / No Precoding');

legend('MU-MIMO, None, TIN',...
       'MU-MIMO, None, R-ML',...
       'Location','southeast');

xlim([0 9]);
ylim([0 75]);

set(gca,'FontSize',12);
