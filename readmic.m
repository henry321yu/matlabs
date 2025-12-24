clc; clear all;close all;

%% ================= 使用者設定 =================
file = "C:\Blastware 10\Event\ascii\ASCII.TXT";   % ← 改成你的實際檔名
Fs_default = 1024;    % 若 metadata 讀不到時的備援 sample rate
%% =============================================

fid = fopen(file,'r');
if fid == -1
    error("無法開啟檔案，請確認路徑與檔名");
end

%% 初始化
meta = struct();      % metadata
MicL1 = [];           % 波形資料
readData = false;     % 是否開始讀取數值

%% 逐行讀檔
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if isempty(line)
        continue;
    end

    % ---------- Metadata ----------
    if startsWith(line,'"')
        line = erase(line,'"');                 % 去掉雙引號
        parts = regexp(line,'^(.*?)\s*:\s*(.*)$','tokens','once');
        if numel(parts) == 2
            key = matlab.lang.makeValidName(strtrim(parts{1}));
            val = strtrim(parts{2});
            meta.(key) = val;
        end

    % ---------- 欄位名稱 ----------
    elseif strcmp(line,'MicL1')
        readData = true;

    % ---------- 波形資料 ----------
    elseif readData
        v = str2double(line);
        if ~isnan(v)
            MicL1(end+1,1) = v; %#ok<SAGROW>
        end
    end
end

fclose(fid);

%% Sample Rate
if isfield(meta,'Sample_Rate')
    Fs = str2double(erase(meta.Sample_Rate,'sps'));
else
    Fs = Fs_default;
end

%% Pre-trigger
preT = 0;
if isfield(meta,'Pre_trigger_Length')
    preT = str2double(erase(meta.Pre_trigger_Length,'sec'));
end

%% 建立時間軸
N = length(MicL1);
t = (0:N-1).' / Fs + preT;

%% 最大值
[maxMic, idx] = max(abs(MicL1));

%% ================ 繪圖 =================
figure;
plot(t, MicL1,'LineWidth',1.2); hold on;
plot(t(idx), MicL1(idx),'ro','LineWidth',1.5);
grid on;

xlabel('Time (sec)');
ylabel('Mic Pressure (Pa)');
title('MiniMate Plus – MicL1 Full Waveform');

text(t(idx), MicL1(idx), ...
    sprintf('  Peak = %.3f Pa @ %.4f s', MicL1(idx), t(idx)), ...
    'VerticalAlignment','bottom');
%% ======================================

%% 印出重點 Metadata
fprintf('\n===== Event Metadata =====\n');
fprintf('Event Type   : %s\n', meta.EventType);
fprintf('Event Date   : %s\n', meta.EventDate);
fprintf('Event Time   : %s\n', meta.EventTime);
fprintf('Trigger      : %s\n', meta.Trigger);
fprintf('Sample Rate  : %.0f sps\n', Fs);
fprintf('Record Time  : %s\n', meta.RecordTime);
fprintf('Pre-trigger  : %.3f sec\n', preT);
fprintf('Mic Peak     : %.3f Pa @ %.4f sec\n', MicL1(idx), t(idx));
fprintf('===========================\n');
