clc;
clear;
close all;

% ====== 使用者設定 ======
samples = 5120;  % 每幀取樣點數（依你的 Acquisition 設定）
% filePath = fullfile('E:\RtDAS\4_DasTcpClient', '2025-10-28_11-37-15.bin');
filePath = fullfile('C:\Users\sgrc - 325\Desktop\git\das\get\das_data', 'raw_data_2025-10-31_14-32-57.bin');
% filePath = fullfile('C:\Users\sgrc - 325\Desktop\git\das\get\das_data', '2025-10-31_14-59-01.bin');
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

% ====== 逐幀讀取並顯示 ======
figure;
for i = 1:frameCount
    % 讀取一幀
    data_raw = fread(fileID, [samples, 1], 'int16=>double');
    
    % 畫圖
    plot(data_raw);
    xlabel('Sampling Points');
    ylabel('ADC Value');
    title(sprintf('Frame %d / %d', i, frameCount));
    % ylim([-8192 8192]);
    grid on;
    
    drawnow; % 即時更新圖形
end

fclose(fileID);
fprintf('✅ 讀取完成，共 %d 幀。\n', frameCount);
