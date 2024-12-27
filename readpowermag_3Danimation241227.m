clear all;close all;

basePath = 'C:\Users\sgrc-325\Desktop\git\pipes\data\p test\';

% fileNames = {'1225_13P.CSV'};  %squar
% fileNames = {'1226_0P.CSV'};  %stop
% fileNames = {'1226_6P.CSV'};  %stop
fileNames = {'1227_0P.CSV'};  %num fix 40 sqaur
% fileNames = {'1227_1P.CSV'};  %stop
% fileNames = {   %multiple data
%     '1224_2P.CSV', ...
%     '1224_3P.CSV', 
%     };
allData = [];

for i = 1:length(fileNames)
    filePath = fullfile(basePath, fileNames{i});
    tempData = readtable(filePath, 'Delimiter', ',', 'ReadVariableNames', false);
    allData = [allData; tempData]; % 將每個表格資料合併
end
%%

moData=allData;

moData=allData(2000:end,:); %1227_0P.CSV
% moData=allData(100:end,:); %1227_1P
% moData=allData(2500:3200,:);
mo=[moData.Var1 moData.Var1 moData.Var1 moData.Var4 moData.Var5 moData.Var6 moData.Var7 moData.Var8 moData.Var9 moData.Var10 moData.Var11 moData.Var12 moData.Var13 moData.Var14 moData.Var15 moData.Var16 moData.Var17 moData.Var18 moData.Var19 moData.Var20];

[No,~]=size(mo);

t= moData.Var1;
time_rtc=moData.Var2;
time_gps=moData.Var3;
m = mo(:,4:6);
acc1 = mo(:,7:9);
gyr1=mo(:,10:12);
lat=mo(:,13);
lon=mo(:,14);
alt=mo(:,15);
sats=mo(:,16);
fix=mo(:,17);
hacc=mo(:,18);
vacc=mo(:,19);
tem=mo(:,20);

ax1=acc1(:,1);ay1=acc1(:,2);az1=acc1(:,3);
gx1=gyr1(:,1);gy1=gyr1(:,2);gz1=gyr1(:,3);
mx=m(:,1);my=m(:,2);mz=m(:,3);
aG=sqrt(ax1.^2+ay1.^2+az1.^2);
GG=sqrt(gyr1(:,1).^2+gyr1(:,2).^2+gyr1(:,3).^2);
mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);
stdd=std(mo)'; %看標準差

% tt = time_rtc;
tt = time_gps;

%%
close all
% 設定投影參數，使用台灣西部二度分帶的 TWD97 座標系
latlon = [lat, lon]; % 經緯度資料組合為矩陣
twd97 = projcrs(3826); % TWD97 座標系

% 使用 projfwd 進行經緯度到座標的轉換
[x, y] = projfwd(twd97, latlon(:,1), latlon(:,2));

x=x*100; %m to cm
y=y*100;
altcm=alt*100;

x=x-x(1);
y=y-y(1);
altcm=altcm-altcm(1);
x1=[];y1=[];altcm1=[];

for i=1:No
    if(fix(i)>3)
        x1(end+1)=x(i);
        y1(end+1)=y(i);
        altcm1(end+1)=altcm(i);
    end
end

fprintf(['標準差 :\n' ...
    'x: %.3fcm,y :%.3fcm,alt: %.3fcm\n' ...
    '測試時長: %.2f 秒\n'], ...
    std(x),std(y),std(altcm),t(end));

fprintf('%s to %s\n', datestr(time_rtc(1)), datestr(time_rtc(end)));

figure(2)
plot(tt,x,'.')
ylabel('TM2 x(cm)');
figure(3)
plot(tt,y,'.')
ylabel('TM2 y(cm)');
figure(4)
plot(tt,altcm1,'.')
ylabel('Altitude(cm)');

figure(5)
plot(x,y,'.')
xlabel('TM2 x(cm)');
ylabel('TM2 y(cm)');
grid on;
axis equal;

figure(6)
plot3(x,y,altcm1,'.')
xlabel('TM2 x(cm)');
ylabel('TM2 y(cm)');
zlabel('Altitude(cm)');
grid on;
axis equal;

% 設定旋轉角度（以度為單位）
theta_deg = 20;  % 旋轉角度設定為 30 度
theta = deg2rad(theta_deg);  % 轉換為弧度
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
rotated_coords = R * [x'; y'];
x_r = rotated_coords(1, :);
y_r = rotated_coords(2, :);

figure(7)
plot(x_r,y_r,'.')
xlabel('TM2 xr(cm)');
ylabel('TM2 yr(cm)');
grid on;
axis equal;

figure(8)
plot3(x_r,y_r,altcm1,'.')
xlabel('TM2 xr(cm)');
ylabel('TM2 yr(cm)');
zlabel('Altitude(cm)');
grid on;
axis equal;

% 初始化圖形
figure(9)
h = plot3(NaN, NaN, NaN, '.');  % 初始化一個空的點圖對象

% 設置圖形屬性
axis tight;  % 自動調整軸範圍
xlabel('TM2 xr(cm)');
ylabel('TM2 yr(cm)');
zlabel('Altitude(cm)');
title('3D Animation');
grid on;
axis equal;

% 循環繪製每一幀
for i = 1:length(x_r)
    % 更新數據
    set(h, 'XData', x_r(1:i), 'YData', y_r(1:i), 'ZData', altcm1(1:i));
%     view([0.0 90.0]) %俯視角
    if mod(i, 1) == 0  % 每50幀更新一次
        drawnow;
    end
end