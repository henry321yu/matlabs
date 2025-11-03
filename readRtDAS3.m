clc;
clear;
close all;

% ====== 使用者設定 ======
samples = 5120;  % 每幀取樣點數
filePath = fullfile('C:\Users\sgrc - 325\Desktop\git\das\get\das_data', ...
                    '2025-10-31_15-48-36.bin');
smoothWindow = 100;  % 平滑視窗長度
% ========================

% 確認檔案存在
if exist(filePath, 'file') ~= 2
    error('❌ 找不到檔案：%s', filePath);
end

% 開啟檔案
fileID = fopen(filePath, 'rb');
if fileID == -1
    error('⚠️ 無法開啟檔案：%s', filePath);
end

% 取得檔案大小（bytes）
fseek(fileID, 0, 'eof');
fileSize = ftell(fileID);
frewind(fileID);

% 計算共有幾幀資料
bytesPerSample = 2;  % int16 = 2 bytes
frameCount = floor(fileSize / (samples * bytesPerSample));

fprintf('📄 檔案大小 = %.2f MB\n', fileSize/1e6);
fprintf('🧩 每幀 %d 點，共 %d 幀\n', samples, frameCount);

% ====== 建立圖形視窗 ======
figure('Name','DAS Data Viewer','NumberTitle','off');
tiledlayout(3,1);

ax1 = nexttile; % 時域
ax2 = nexttile; % FFT 震幅
ax3 = nexttile; % FFT 相位

for i = 1:frameCount
    % === 讀取一幀 ===
    data_raw = fread(fileID, [samples, 1], 'int16=>double');
    if isempty(data_raw)
        break;
    end

    % === FFT 分析 ===
    data_fft = fft(data_raw);
    amplitude = abs(data_fft);
    phase = angle(data_fft);

    % === 平滑處理 ===
    amplitude_smooth = smoothdata(amplitude, 'movmean', smoothWindow);
    phase_smooth = smoothdata(phase, 'movmean', smoothWindow);

    % === 繪圖 ===
    % 時域
    axes(ax1);
    plot(data_raw, 'b');
    xlabel('Sampling Points');
    ylabel('ADC Value');
    title(sprintf('Frame %d / %d (時域訊號)', i, frameCount));
    grid on;

    % FFT 震幅
    axes(ax2);
    plot(amplitude_smooth, 'r');
    xlabel('Frequency Bin');
    ylabel('Amplitude');
    title(sprintf('Frame %d / %d (FFT 震幅, smooth=%d)', i, frameCount, smoothWindow));
    grid on;

    % FFT 相位
    axes(ax3);
    plot(phase_smooth, 'k');
    xlabel('Frequency Bin');
    ylabel('Phase (radians)');
    title(sprintf('Frame %d / %d (FFT 相位, smooth=%d)', i, frameCount, smoothWindow));
    grid on;

    drawnow; % 即時更新圖形
end

fclose(fileID);
fprintf('✅ 所有幀播放完成，共 %d 幀。\n', frameCount);
