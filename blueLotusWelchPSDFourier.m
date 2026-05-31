% Generated on: 31-May-2026 14:05:06
data = OpenBCIRAW20231203220519;
data = data(7:end, :);
data = table2array(data);
data = data(3:end, :);

lotEEG = data;            % 16 EEG channels (cols 2..17)

Fs = 125;                          % sampling rate (Hz)

%% Welch PSD parameters
windowLen = 3285;                  % ~26 s window
overlap   = 1642;                  % 50% overlap
nfft      = 3285;
winFunc   = hamming(windowLen, 'periodic');

%% Compute PSD for all 16 channels at once
nChan = size(lotEEG, 2);
[P, F] = pwelch(lotEEG(:,1), winFunc, overlap, nfft, Fs);  % first call -> F
P_dB        = zeros(numel(F), nChan);
P_dB(:,1)   = pow2db(P);
for c = 2:nChan
    P = pwelch(lotEEG(:,c), winFunc, overlap, nfft, Fs);
    P_dB(:,c) = pow2db(P);
end

%% Plot
chanNames = {'Fp1','Fp2','C3','C4','T5','T6','O1','O2', ...
             'F7','F8','F3','F4','T3','T4','P3','P4'};

fig = figure('Color','white','Position',[60 60 1400 700]);
plot(F, P_dB, 'LineWidth', 0.9);
xlabel('Frequency (Hz)','FontSize',12);
ylabel('Power Spectrum (dB)','FontSize',12);
title('Blue Lotus EEG -- Welch PSD, 16 channels','FontSize',13);
xlim([0 62.5]);                      % cuts the 60 Hz line-noise tail
grid on;
legend(chanNames, 'Location','eastoutside','FontSize',9);

%% Save (drop into your book's images folder)
exportgraphics(fig, 'lotusFourierFreq.jpg', 'Resolution', 200);
