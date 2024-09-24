% 定義資料夾路徑
folder_path = 'C:\Users\sgrc-325\Desktop\py\rasp log';

% 初始化一個空的 table 來儲存所有的資料
all_data = table();

% 用來儲存顏色的字典
color_map = containers.Map('KeyType', 'char', 'ValueType', 'any');

% 獲取資料夾中所有.txt檔案
files = dir(fullfile(folder_path, '*.txt'));

% 遍歷所有檔案
for i = 1:length(files)
    filename = files(i).name;
    file_path = fullfile(folder_path, filename);
    
    % 提取日期
    date_str = extractBetween(filename, '_', '_');
    date = datetime(date_str, 'InputFormat', 'yyyyMMdd');

    % 隨機生成顏色
    if ~isKey(color_map, datestr(date, 'yyyy-mm-dd'))
        color_map(datestr(date, 'yyyy-mm-dd')) = rand(1, 3); % 隨機顏色
    end
    color = color_map(datestr(date, 'yyyy-mm-dd'));

    % 初始化一個暫存表來存儲數據
    temp_data = [];

    % 打開檔案逐行讀取
    fid = fopen(file_path, 'r');
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, 'END_OF_FILE')
            break;
        end

        split_line = strsplit(strtrim(line));
        
        % 檢查是否只有9個欄位，若是則補上'volt'和'current'為0
        if length(split_line) == 9
            split_line = [split_line, {'0', '0'}]; % 增加 'volt' 和 'current' 欄位為 0
        end
        
        temp_data = [temp_data; split_line];
    end
    fclose(fid);

    % 如果有資料，則將其轉換為 table
    if ~isempty(temp_data)
        data = cell2table(temp_data);
        data.Properties.VariableNames = {'UTC+8', 'ax', 'ay', 'az', 'lat', 'lon', 'altitude', 'gps_mode', 'temperature', 'volt', 'current'};
        
        % 將時間欄轉換為datetime格式
        data.UTC8 = datetime(data.UTC8, 'InputFormat', 'HH:mm:ss.SSS', 'Format', 'HH:mm:ss.SSS');

        % 丟棄無法轉換時間的行
        data(isnan(data.UTC8), :) = [];

        % 將所有除了'UTC+8'的欄位轉換為數字
        for j = 2:width(data)
            data{:, j} = str2double(data{:, j});
        end

        % 丟棄任何包含 NaN 的行
        data(any(isnan(data{:, 2:end}), 2), :) = [];

        % 將 lat 和 lon 轉換為 TWD97 (EPSG:3826)
        [twd97_x, twd97_y] = projfwd(projcrs(3826), data.lat, data.lon);
        data.twd97_x = twd97_x;
        data.twd97_y = twd97_y;

        % 將日期新增到資料中
        data.date = repmat(date, height(data), 1);

        % 將資料彙整起來
        all_data = [all_data; data];
    end
end

% 確保數據按時間順序排序
all_data = sortrows(all_data, 'UTC_8');

% 繪製圖表，根據日期使用不同顏色
figure; hold on;
for date_key = keys(color_map)
    date_data = all_data(all_data.date == datetime(date_key{1}), :);
    plot(date_data.UTC_8, date_data.volt, '.', 'DisplayName', date_key{1}, 'Color', color_map(date_key{1}));
end
ylim([5, 20]);
datetick('x', 'HH:MM:SS', 'keepticks');
xlabel('Time (UTC+8)');
ylabel('Volt');
title('Volt vs Time');
legend('show');
hold off;

% 繪製 current 圖表
figure; hold on;
for date_key = keys(color_map)
    date_data = all_data(all_data.date == datetime(date_key{1}), :);
    plot(date_data.UTC_8, date_data.current, '.', 'DisplayName', date_key{1}, 'Color', color_map(date_key{1}));
end
datetick('x', 'HH:MM:SS', 'keepticks');
xlabel('Time (UTC+8)');
ylabel('Current');
title('Current vs Time');
legend('show');
hold off;

% 繪製 gps_mode 圖表
figure; hold on;
for date_key = keys(color_map)
    date_data = all_data(all_data.date == datetime(date_key{1}), :);
    plot(date_data.UTC_8, date_data.gps_mode, '.', 'DisplayName', date_key{1}, 'Color', color_map(date_key{1}));
end
datetick('x', 'HH:MM:SS', 'keepticks');
xlabel('Time (UTC+8)');
ylabel('GPS Mode');
title('GPS Mode vs Time');
legend('show');
hold off;

% 繪製 temperature 圖表
figure; hold on;
for date_key = keys(color_map)
    date_data = all_data(all_data.date == datetime(date_key{1}), :);
    plot(date_data.UTC_8, date_data.temperature, '.', 'DisplayName', date_key{1}, 'Color', color_map(date_key{1}));
end
datetick('x', 'HH:MM:SS', 'keepticks');
xlabel('Time (UTC+8)');
ylabel('Temperature');
title('Temperature vs Time');
legend('show');
hold off;

% 只選擇 gps_mode 為 4 的資料
gps_mode_4_data = all_data(all_data.gps_mode == 4, :);

% 繪製 twd97_x vs UTC+8 (僅限 gps_mode 為 4 的數據)
figure; hold on;
for date_key = keys(color_map)
    date_gps_data = gps_mode_4_data(gps_mode_4_data.date == datetime(date_key{1}), :);
    plot(date_gps_data.UTC_8, date_gps_data.twd97_x, '.', 'DisplayName', date_key{1}, 'Color', color_map(date_key{1}));
end
datetick('x', 'HH:MM:SS', 'keepticks');
xlabel('Time (UTC+8)');
ylabel('TWD97 X');
title('TWD97 X vs Time (GPS Mode=4)');
legend('show');
hold off;

% 繪製 twd97_y vs UTC+8 (僅限 gps_mode 為 4 的數據)
figure; hold on;
for date_key = keys(color_map)
    date_gps_data = gps_mode_4_data(gps_mode_4_data.date == datetime(date_key{1}), :);
    plot(date_gps_data.UTC_8, date_gps_data.twd97_y, '.', 'DisplayName', date_key{1}, 'Color', color_map(date_key{1}));
end
datetick('x', 'HH:MM:SS', 'keepticks');
xlabel('Time (UTC+8)');
ylabel('TWD97 Y');
title('TWD97 Y vs Time (GPS Mode=4)');
legend('show');
hold off;

% 檢查結果
disp(all_data(:, {'UTC_8', 'lat', 'lon', 'twd97_x', 'twd97_y', 'volt', 'current'}));
