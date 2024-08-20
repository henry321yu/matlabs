clear all

% 輸入編碼器量測結果inputEncoder，
% 解算出 速度歷時turning_velocity 與 累積位移量歷時Encoder_X

% inputEncoder = load('inputEncoder.txt'); % 編碼器量測數據 須設定
mo = load('C:\Users\sgrc-325\OneDrive - 國立成功大學 National Cheng Kung University\桌面\temp\log526.csv');
inputEncoder = [mo(:,1) mo(:,21)];
Encoder_threshold=30; % 參數須設定
circumference=0.66; % 輪子周長(公尺) 須設定

[N a2]=size(inputEncoder);

peak_max=inputEncoder(1,2);
peak_min=inputEncoder(1,2);

condition=0;
NO_peak_max=0;
NO_peak_min=0;
peak_max_i=zeros(N,1);
peak_min_i=zeros(N,1);
for i=1:N
    
    if condition==0 % 剛開始的時候
        if peak_max<inputEncoder(i,2)
            peak_max=inputEncoder(i,2);
            peak_max_i_temp=i;
        end
        if peak_min>inputEncoder(i,2)
            peak_min=inputEncoder(i,2);
            peak_min_i_temp=i;
        end
        if peak_max-peak_min>Encoder_threshold......
                && inputEncoder(i,2)<inputEncoder(1,2)
            condition=2;
            NO_peak_max=NO_peak_max+1;
            peak_max_i(NO_peak_max,1)=peak_max_i_temp;
            peak_max=peak_min;
        end
        if abs(inputEncoder(i,2)-peak_min)>Encoder_threshold......
                && inputEncoder(i,2)>inputEncoder(1,2)
            condition=1;
            NO_peak_min=NO_peak_min+1;
            peak_min_i(NO_peak_min,1)=peak_min_i_temp;
            peak_min=peak_max;
        end
    end
    
    if condition==1 % 編碼器的磁力強度在平均值之上 正在找尋波峰值
        if peak_max<inputEncoder(i,2)
            peak_max=inputEncoder(i,2);
            peak_max_i_temp=i;
        end
        if peak_min>inputEncoder(i,2)
            peak_min=inputEncoder(i,2);
            peak_min_i_temp=i;
        end
        if peak_max-peak_min>Encoder_threshold
            condition=2;
            NO_peak_max=NO_peak_max+1;
            peak_max_i(NO_peak_max,1)=peak_max_i_temp;
            peak_max=peak_min;
        end
    end
    
    if condition==2 % 編碼器的磁力強度在平均值之下 正在找尋波谷值
        if peak_max<inputEncoder(i,2)
            peak_max=inputEncoder(i,2);
            peak_max_i_temp=i;
        end
        if peak_min>inputEncoder(i,2)
            peak_min=inputEncoder(i,2);
            peak_min_i_temp=i;
        end
        if peak_max-peak_min>Encoder_threshold
            condition=1;
            NO_peak_min=NO_peak_min+1;
            peak_min_i(NO_peak_min,1)=peak_min_i_temp;
            peak_min=peak_max;
        end
    end
    
end
peak_max_i=peak_max_i(1:NO_peak_max,1);
peak_min_i=peak_min_i(1:NO_peak_min,1);

turning_flag=zeros(N,1);
turning_degree=zeros(N,1);
turning_velocity=zeros(N,1);
dis_degree=zeros(N,1);
Encoder_X=zeros(N,1);

turning_flag(peak_max_i(:,1),1)=1;
turning_flag(peak_min_i(:,1),1)=-1;

turning_condition=0;
for i=1:N % 磁力由正轉負
    if turning_flag(i,1)==1
        turning_condition=1;
        temp1=i;
    end
    if turning_condition==1 && turning_flag(i,1)==-1
        temp2=i;
        turning_degree(temp1:temp2,1)=.....
            (inputEncoder(temp1:temp2,2)-inputEncoder(temp2,2))*180.....
            /(inputEncoder(temp1,2)-inputEncoder(temp2,2))-90;
        turning_flag(temp1+1:temp2-1,1)=0.5;
    end
end
turning_condition=0;
for i=1:N % 磁力由負轉正
    if turning_flag(i,1)==-1
        turning_condition=1;
        temp1=i;
    end
    if turning_condition==1 && turning_flag(i,1)==1
        temp2=i;
        turning_degree(temp1:temp2,1)=.....
            (inputEncoder(temp1:temp2,2)-inputEncoder(temp1,2))*180.....
            /(inputEncoder(temp2,2)-inputEncoder(temp1,2))-90;
        turning_flag(temp1+1:temp2-1,1)=-0.5;
    end
end

for i=1:N-1
    if  turning_flag(i,1)==-0.5 || turning_flag(i+1,1)==-0.5 % 磁力由負轉正
        dis_degree(i,1)=turning_degree(i+1,1)-turning_degree(i,1); % 角度差
    end
    if  turning_flag(i,1)==+0.5 || turning_flag(i+1,1)==+0.5 % 磁力由正轉負
        dis_degree(i,1)=turning_degree(i,1)-turning_degree(i+1,1); % 角度差
    end
end

for i=2:N-1
    ave_degree=(dis_degree(i-1,1)+dis_degree(i,1))/2; % 平均角度差
    dis_time=(inputEncoder(i+1,1)-inputEncoder(i-1,1))/2; % 時間差
    turning_velocity(i,1)=ave_degree/360*circumference/dis_time; % 速度歷時
    Encoder_X(i+1,1)=Encoder_X(i,1)+ave_degree/360*circumference; % 累積位移量歷時
end

figure(1)
plot(inputEncoder(:,1),inputEncoder(:,2));
hold on
plot(inputEncoder(peak_max_i(:,1),1),inputEncoder(peak_max_i(:,1),2),'*');
hold on
plot(inputEncoder(peak_min_i(:,1),1),inputEncoder(peak_min_i(:,1),2),'+');
figure(2)
plot(inputEncoder(:,1),turning_degree(:,1));
figure(3)
plot(inputEncoder(:,1),turning_velocity(:,1));
figure(4)
plot(inputEncoder(:,1),Encoder_X(:,1));
