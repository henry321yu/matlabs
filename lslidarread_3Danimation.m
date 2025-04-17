clear all;
% lidar=load(['C:\Users\sgrc-325\Desktop\git\lidar\lslidar.txt']);
% 讀取 CSV 檔案
% data = readtable(['C:\Users\sgrc-325\Desktop\git\lidar\LidarType_C_v5_0_2025-04-15-16-38_1Frame.csv']); %20hz 
% data = readtable(['C:\Users\sgrc-325\Desktop\git\lidar\LidarType_C_v5_0_2025-04-15-16-38_2Frame.csv']);
% data = readtable(['C:\Users\sgrc-325\Desktop\git\lidar\LidarType_C_v5_0_2025-04-15-16-38_200Frame.csv']);

data = readtable(['C:\Users\sgrc-325\Desktop\git\lidar\LidarType_C_v5_0_2025-04-16-11-58_1Frame.csv']); %5hz
% data = readtable(['C:\Users\sgrc-325\Desktop\git\lidar\LidarType_C_v5_0_2025-04-16-11-58_100Frame.csv']);

close all;
% 取出座標欄位
x = data.Points_X;
y = data.Points_Y;
z = data.Points_Z;
ang=data.H_angle;
dis=data.Distance1;
int=data.Intensity1;
ts=data.Timestamp_s;
tns=data.Timestamp_ns;

t=(ts-ts(1))+(tns-tns(1))*1e-9;
[N,~]=size(t);
f=N/t(end)/1000;
fprintf('time length: %.3f s\n',t(end));
fprintf('f: %.2f kHz\n',f);


% 繪製 3D 點雲圖
figure(1);
plot3(x, y, z, '.');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('LiDAR Point Cloud');
grid on;
axis equal;

view([0 90])

figure(2)
plot(t,ang,'.')
xlim([0 0.5])
xlabel('time(s)');
ylabel('angle');


figure(3)
h = plot3(NaN, NaN, NaN, '.');  % 初始化一個空的點圖對象

% 設置圖形屬性
axis tight;  % 自動調整軸範圍
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('LiDAR Point Cloud Animation');
grid on;
axis equal;

% 循環繪製每一幀
for i = 1:length(x)
    % 更新數據
    set(h, 'XData', x(1:i), 'YData', y(1:i), 'ZData', z(1:i));

    % view([-180 60]) %正視仰角
    % view([0 90]) %俯視角
    % view([-20 10]) %斜視角
    if mod(i, 1000) == 0  % 每k幀更新一次
        drawnow;
    end
end
