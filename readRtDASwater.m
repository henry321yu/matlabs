clc
clear
close all
samples = 640;
% filePath = fullfile('C:\Users\sgrc - 325\Desktop\git\das\get\das_data', ...
%                     '2025-10-31_16-33-02.bin');

% filePath = fullfile('C:\Users\sgrc - 325\Desktop\git\das\get\das_data', ...
%                     'raw_data_2025-10-31_16-32-47.bin');



% offical demo
% filePath = fullfile('C:\Users\sgrc - 325\Desktop\temp_sync', ...
%                     'raw_data_2025-11-03_14-59-17.bin');
filePath = fullfile('C:\Users\sgrc - 325\Desktop\temp_sync', ...
                    'demod_data_2025-11-03_15-02-33.bin');

fileID = fopen(filePath, 'rb');
phase = fread(fileID,[samples, inf],'float=>float');
fclose(fileID);
imagesc(phase',[-0.5,0.5])
colormap('jet')
xlabel('location(point)');
ylabel('pulse(n)');
title('Das-Pro phase waterfall')