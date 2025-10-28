clc;
clear;
close all;

% 依據你的取樣設定輸入樣本數 (samples)
samples = 2560;  % ← 需與當時 Acquisition Parameter 中設定的 sampling length 一致

% 指定檔案路徑
filePath = 'C:\Users\sgrc-325\Desktop\2025-10-28_11-37-15.bin';

% 開啟檔案
fileID = fopen(filePath, 'rb');
if fileID == -1
    error('無法開啟檔案：%s', filePath);
end

% 讀取資料 (int16格式轉為double)
data_raw = fread(fileID, [samples, inf], 'int16=>double');
fclose(fileID);

% 顯示資料資訊
fprintf('成功讀取資料，大小：%d x %d\n', size(data_raw,1), size(data_raw,2));

% 顯示第一筆波形
figure;
plot(data_raw(:,1));
xlabel('Sampling Points');
ylabel('ADC Value');
title('Raw Data Waveform');
grid on;
