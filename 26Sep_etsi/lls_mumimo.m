clear; clc; rng(0);

%% Parameters
snrVec = 0:2:9;
Nt = 4;  Nr = 2;
nFrames = 10;
carrier = nrCarrierConfig;
carrier.NSizeGrid = 51;
carrier.SubcarrierSpacing = 30;
slotsPerFrame = carrier.SlotsPerFrame;      % 20 at 30 kHz
Tslot = 1e-3 / (carrier.SubcarrierSpacing/15);   % 0.5 ms

%% PDSCH (two co-scheduled UEs, orthogonal DM-RS ports)
pdschUE1 = nrPDSCHConfig;
pdschUE1.PRBSet = 0:50;  pdschUE1.SymbolAllocation = [2 12];
pdschUE1.MappingType = "A";  pdschUE1.Modulation = "QPSK";
pdschUE1.NumLayers = 1;  pdschUE1.RNTI = 101;
pdschUE1.DMRS.DMRSPortSet = 0;

pdschUE2 = pdschUE1;
pdschUE2.RNTI = 102;
pdschUE2.DMRS.DMRSPortSet = 2;

%% OFDM / channel
ofdmInfo = nrOFDMInfo(carrier);
Fs = ofdmInfo.SampleRate;  Nfft = ofdmInfo.Nfft;

tdl = nrTDLChannel;
tdl.DelayProfile = 'TDL-C';
tdl.DelaySpread = 300e-9;
tdl.NumTransmitAntennas = Nt;
tdl.NumReceiveAntennas = Nr;
tdl.SampleRate = Fs;
tdl.ChannelResponseOutput = "ofdm-response";
chInfo = info(tdl);
maxChDelay = ceil(max(chInfo.PathDelays*Fs)) + chInfo.ChannelFilterDelay;

% QPSK constellation in bit order [00 01 10 11] (for ML detector)
bt = [0 0; 0 1; 1 0; 1 1];
const = nrSymbolModulate(reshape(bt.',[],1),"QPSK");
slice = @(x) reshape([real(x(:)).'<0; imag(x(:)).'<0], [], 1);   % TS 38.211 QPSK

%% Results
thrTIN = zeros(size(snrVec)); thrML = thrTIN; thrSU = thrTIN;
berTIN = thrTIN; berML = thrTIN; berSU = thrTIN;

for isnr = 1:numel(snrVec)
    snr = snrVec(isnr);
    N0  = 10^(-snr/10);                       % RE-level noise variance
    errTIN = 0; errML = 0; errSU = 0; nBits = 0; nSlots = 0;
    reset(tdl);

    for iframe = 1:nFrames
        for islot = 1:slotsPerFrame
            carrier.NSlot = islot-1;          % 0..19

            %% Tx grid
            txGrid = nrResourceGrid(carrier,Nt);
            nRE = size(txGrid,1)*size(txGrid,2);

            [ind1,info1] = nrPDSCHIndices(carrier,pdschUE1);
            [ind2,info2] = nrPDSCHIndices(carrier,pdschUE2);
            assert(isequal(ind1,ind2),'UEs must share data REs');

            bits1 = randi([0 1],info1.G,1);
            bits2 = randi([0 1],info2.G,1);
            s1 = nrSymbolModulate(bits1,"QPSK");
            s2 = nrSymbolModulate(bits2,"QPSK");

            txGrid(ind1)       = s1;          % UE1 -> antenna 1
            txGrid(ind2 + nRE) = s2;          % UE2 -> antenna 2

            txGrid(nrPDSCHDMRSIndices(carrier,pdschUE1))       = nrPDSCHDMRS(carrier,pdschUE1);
            txGrid(nrPDSCHDMRSIndices(carrier,pdschUE2) + nRE) = nrPDSCHDMRS(carrier,pdschUE2);

            %% OFDM + channel
            txWave = nrOFDMModulate(carrier,txGrid);
            txWave = [txWave; zeros(maxChDelay,Nt)];
            [rxWave,ofdmResp,tOff] = tdl(txWave,carrier);

            if tOff > 0
                rxWave = rxWave(1+tOff:end,:);
            elseif tOff < 0
                rxWave = [zeros(-tOff,Nr,'like',rxWave); rxWave];
            end

            %% Noise: RE variance N0  =>  time-domain variance N0/Nfft
            noise = sqrt(N0/(2*Nfft))*(randn(size(rxWave))+1i*randn(size(rxWave)));
            rxGrid = nrOFDMDemodulate(carrier,rxWave+noise);

            %% Extract data REs and their channels
            [y,H] = nrExtractResources(ind1,rxGrid,ofdmResp);   % y: Nre x Nr, H: Nre x Nr x Nt
            H1 = H(:,:,1);                    % Nre x Nr, UE1 (desired)
            H2 = H(:,:,2);                    % Nre x Nr, UE2 (interference)

            %% (a) TIN: MMSE-IRC treating UE2 as noise, per RE
            R11 = abs(H2(:,1)).^2 + N0;
            R22 = abs(H2(:,2)).^2 + N0;
            R12 = H2(:,1).*conj(H2(:,2));
            det = R11.*R22 - abs(R12).^2;
            w = [ (R22.*H1(:,1) - R12.*H1(:,2)), ...
                  (R11.*H1(:,2) - conj(R12).*H1(:,1)) ] ./ det;   % Rin \ h1
            xTIN = sum(conj(w).*y,2) ./ sum(conj(w).*H1,2);
            bHatTIN = slice(xTIN);

            %% (b) Joint ML over 16 (s1,s2) pairs, per RE
            best = inf(size(y,1),1); i1 = ones(size(best)); i2 = i1;
            for a = 1:4
                for b = 1:4
                    d = sum(abs(y - H1*const(a) - H2*const(b)).^2,2);
                    u = d < best;
                    best(u) = d(u); i1(u) = a; i2(u) = b;
                end
            end
            bHatML = reshape(bt(i1,:).',[],1);

            % %% (c) Genie interference-free bound (UE2 removed, MRC)
            % ySU = y - H2.*repmat(s2,1,Nr);            % Nre x Nr
            % xSU = sum(conj(H1).*ySU,2) ./ sum(abs(H1).^2,2);
            % bHatSU = slice(xSU);

            %% Bookkeeping
            errTIN = errTIN + sum(bHatTIN ~= bits1);
            errML  = errML  + sum(bHatML  ~= bits1);
            % errSU  = errSU  + sum(bHatSU  ~= bits1);
            nBits  = nBits + numel(bits1);
            nSlots = nSlots + 1;
        end
    end

    T = nSlots*Tslot;
    berTIN(isnr) = errTIN/nBits;  
    thrTIN(isnr) = (nBits-errTIN)/T/1e6;

    berML(isnr)  = errML/nBits;   
    thrML(isnr)  = (nBits-errML)/T/1e6;

    % berSU(isnr)  = errSU/nBits;   thrSU(isnr)  = (nBits-errSU)/T/1e6;
    % fprintf('SNR=%2d dB | BER TIN=%.3e ML=%.3e SU=%.3e | Thr TIN=%6.2f ML=%6.2f SU=%6.2f Mbps\n', ...
    %     snr,berTIN(isnr),berML(isnr),berSU(isnr),thrTIN(isnr),thrML(isnr),thrSU(isnr));

    fprintf('SNR=%2d dB | Thr TIN=%6.2f ML=%6.2f Mbps\n', ...
        snr,thrTIN(isnr),thrML(isnr));
end

%% Plots
% figure;
% semilogy(snrVec,max(berTIN,1e-7),'bo-',snrVec,max(berML,1e-7),'rs-', ...
%          snrVec,max(berSU,1e-7),'k^--','LineWidth',2);
% grid on; xlabel('SNR (dB)'); ylabel('UE1 uncoded BER');
% legend('TIN (MMSE-IRC)','Joint ML','Interference-free','Location','southwest');

figure;
plot(snrVec,thrTIN,'ro-',snrVec,thrML,'bo-','LineWidth',2);
grid on; xlabel('SNR (dB)'); ylabel('Throughput (Mbps)');
legend('MU-MIMO, None, TIN','MU-MIMO, None, Joint ML','northwest');