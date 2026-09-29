clear;
clc;
rng(0);

%% Parameters
snrVec = 0:2:9;

Nt = 4;
Nr = 2;

nFrames = 10; %100
nSlotsPerFrame = 20; %20

carrier = nrCarrierConfig;
carrier.NSizeGrid = 51;
carrier.SubcarrierSpacing = 30;

%% PDSCH
% UE1
pdschUE1 = nrPDSCHConfig;
pdschUE1.PRBSet           = 0:50;
pdschUE1.SymbolAllocation = [2 12];
pdschUE1.MappingType      = "A";
pdschUE1.Modulation       = "QPSK";
pdschUE1.NumLayers        = 1;
pdschUE1.RNTI             = 101;
pdschUE1.DMRS.DMRSPortSet = 0;       % AP1000

% UE2
pdschUE2 = nrPDSCHConfig;
pdschUE2.PRBSet           = 0:50;
pdschUE2.SymbolAllocation = [2 12];
pdschUE2.MappingType      = "A";
pdschUE2.Modulation       = "QPSK";
pdschUE2.NumLayers        = 1;
pdschUE2.RNTI             = 102;
pdschUE2.DMRS.DMRSPortSet = 2;       % AP1002

%% OFDM
ofdmInfo = nrOFDMInfo(carrier);
Fs = ofdmInfo.SampleRate;

%% Channel
tdl = nrTDLChannel;

tdl.DelayProfile = 'TDL-C';
tdl.DelaySpread = 300e-9;

tdl.NumTransmitAntennas = Nt;
tdl.NumReceiveAntennas = Nr;
tdl.SampleRate = Fs;
tdl.ChannelResponseOutput = "ofdm-response";

%% Throughput Results
thrTIN = zeros(size(snrVec));
thrRML = zeros(size(snrVec));

%% Main SNR loop

for isnr = 1:length(snrVec)

    snr = snrVec(isnr);

    totalBitsTIN = 0;
    totalBitsRML = 0;

    totalTime = 0;

    reset(tdl);
    for iframe = 1:nFrames
        for islot = 1:nSlotsPerFrame

            carrier.NSlot = islot;

            %% Resource Grid
            % generates an empty carrier slot resource grid for the specified number of antennas
            txGrid = nrResourceGrid(carrier,Nt);
            
            %% PDSCH Resources
            [pdschInd1,pdschInfo1] = nrPDSCHIndices(carrier,pdschUE1);
            [pdschInd2,pdschInfo2] = nrPDSCHIndices(carrier,pdschUE2);
            
            bitsUE1 = randi([0 1],pdschInfo1.G,1);
            symUE1 = nrSymbolModulate(bitsUE1,"QPSK");
            bitsUE2 = randi([0 1],pdschInfo2.G,1);
            symUE2 = nrSymbolModulate(bitsUE2,"QPSK");

            %% MU-MIMO Mapping
            nRE = numel(txGrid(:,:,1));

            txGrid(pdschInd1) = symUE1;
            txGrid(pdschInd2 + nRE) = symUE2;

            %% DMRS UE1
            % Generate PDSCH DM-RS symbols
            dmrsSym1 = nrPDSCHDMRS(carrier,pdschUE1);
            dmrsInd1 = nrPDSCHDMRSIndices(carrier,pdschUE1);

            txGrid(dmrsInd1) = dmrsSym1;

            %% DMRS UE2
            dmrsSym2 = nrPDSCHDMRS(carrier,pdschUE2);
            dmrsInd2 = nrPDSCHDMRSIndices(carrier,pdschUE2);

            txGrid(dmrsInd2 + numel(txGrid(:,:,1))) = dmrsSym2;

            %% OFDM
            txWaveform = nrOFDMModulate(carrier,txGrid);

            %% Channel
            %[rxWaveform,pathGains] = tdl(txWaveform);
            [rxWaveform,ofdmResponse,timingOffset] = tdl(txWaveform,carrier);

            %% Apply timing synchronization
            if timingOffset > 0
                rxWaveform = rxWaveform(1 + timingOffset:end, :);
            elseif timingOffset < 0
                rxWaveform = [zeros(-timingOffset,size(rxWaveform,2),'like',rxWaveform); rxWaveform];
            end

            %pathFilters = getPathFilters(tdl);

            %% Perfect CSIR and CSIT
            % Hperfect = nrPerfectChannelEstimate(carrier, ...
            %     pathGains, pathFilters);

            %% Noise
            N0 = 10^(-snr/10);

            noise = sqrt(N0/2)*(randn(size(rxWaveform)) + 1i*randn(size(rxWaveform)));
            rxWaveform = rxWaveform + noise;

            %% OFDM Demodulation
            rxGrid = nrOFDMDemodulate(carrier,rxWaveform);

            %% Extract PDSCH REs

            [pdschRx,pdschH] = nrExtractResources( ...
                pdschInd1,...
                rxGrid,...
                ofdmResponse); % Hperfect

            %% Average channel
            Havg = squeeze(mean(pdschH,1));

            %% Desired user channel
            H1 = Havg(:,1);

            %% Interfering user channel
            H2 = Havg(:,2);

            %% -------------------------------------------------
            %% TIN Detector
            %% -------------------------------------------------

            Irx = eye(Nr);

            % MMSE/TIN combiner for UE1, treating UE2 as interference
            Wtin = ((H1 * H1' + H2 * H2' + N0 * Irx) \ H1);   % Nr x 1

            % Equalized UE1 symbols
            xTIN = pdschRx * conj(Wtin);            % Nre x 1

            % QPSK hard-decision demodulation
            bitsHatTIN = qamdemod(xTIN,4, ...
                'gray', ...
                'UnitAveragePower',true, ...
                'OutputType','bit');
            bitsHatTIN = bitsHatTIN(:);

            % Compare against transmitted UE1 bits
            refBitsTIN = bitsUE1(1:numel(bitsHatTIN));
            nCorrectTIN = sum(bitsHatTIN == refBitsTIN);

            % Throughput as correctly recovered bits
            totalBitsTIN = totalBitsTIN + nCorrectTIN;


            % Wtin = (H1'*H1 + N0)\H1';
            % xTIN = Wtin*pdschRx';
            % 
            % bitsHatTIN = real(xTIN)<0;
            % 
            % nCorrectTIN = sum(bitsHatTIN == bitsUE1(1:length(bitsHatTIN)));
            % 
            % if nCorrectTIN > 0.99*length(bitsHatTIN)
            %     totalBitsTIN = totalBitsTIN + length(bitsUE1);
            % end

            totalTime = totalTime + 0.5e-3;
            % 0.5e-3 is the duration of one NR slot when SubcarrierSpacing = 30 kHz.
            % So this line is counting simulated time per slot.
            
        end
    end

    thrTIN(isnr) = totalBitsTIN/totalTime/1e6;

    fprintf('SNR=%2d dB TIN=%6.2f Mbps\n', snr, thrTIN(isnr));
end

%% Plot
figure;

plot(snrVec,thrTIN,'bo-','LineWidth',2);
grid on;
xlabel('SNR (dB)');
ylabel('Average Throughput (Mbps)');
legend('TIN','Location','northwest');
