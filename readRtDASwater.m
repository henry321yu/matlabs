clc
clear
close all
samples = 640;
filePath = fullfile('C:\Users\sgrc - 325\Desktop\git\das\get\das_data', ...
                    '2025-10-31_16-33-02.bin');

% filePath = fullfile('C:\Users\sgrc - 325\Desktop\git\das\get\das_data', ...
%                     'raw_data_2025-10-31_16-32-47.bin');

fileID = fopen(filePath, 'rb');
phase = fread(fileID,[samples, inf],'float=>float');
fclose(fileID);
imagesc(phase',[-0.5,0.5])
colormap('jet')
xlabel('location(point)');
ylabel('pulse(n)');
title('Das-Pro phase waterfall')