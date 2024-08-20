clear; close all;


fitdata=xlsread('fit data\9 30 fit data.xlsx');
Ba=fitdata(1:3,1);Ta=fitdata(1:3,2:4);Ka=fitdata(1:3,5:7);
Bg=fitdata(5:7,1);Tg=fitdata(5:7,2:4);Kg=fitdata(5:7,5:7);
Bm=fitdata(9:11,1);Tm2a=fitdata(9:11,2:4);

% 時間  加速度(m/s/s)XYZ 角速度XYZ 磁力XYZ 放在data.txt 中
data=load('fit data\9 25 new fit code.csv');

acc=data(:,2:4);
gyr=data(:,5:7);
mag=data(:,8:10);
Aa=sqrt(acc(:,1).^2+acc(:,2).^2+acc(:,3).^2);
Ag=sqrt(gyr(:,1).^2+gyr(:,2).^2+gyr(:,3).^2);
Am=sqrt(mag(:,1).^2+mag(:,2).^2+mag(:,3).^2);


% a=1/2048;g=1/16.4;m=0.6;
a=1;g=1;m=1;
sm=25;
data(:,2)=smooth(data(:,2),sm)*a; %acc
data(:,3)=smooth(data(:,3),sm)*a;
data(:,4)=smooth(data(:,4),sm)*a;
data(:,5)=smooth(data(:,5),sm)*g; %gyr
data(:,6)=smooth(data(:,6),sm)*g;
data(:,7)=smooth(data(:,7),sm)*g;
data(:,8)=smooth(data(:,8),sm)*m; %mag
data(:,9)=smooth(data(:,9),sm)*m;
data(:,10)=smooth(data(:,10),sm)*m;
threshold=2.7; %3.2(2 possible) 5.6(1 possible) min 2.7   9 30 box test 

% acc2=data(:,2:4);
% gyr2=data(:,5:7);
% mag2=data(:,8:10);
% Aa2=sqrt(acc2(:,1).^2+acc2(:,2).^2+acc2(:,3).^2);
% Ag2=sqrt(gyr2(:,1).^2+gyr2(:,2).^2+gyr2(:,3).^2);
% Am2=sqrt(mag2(:,1).^2+mag2(:,2).^2+mag2(:,3).^2);

% [Ta,Ka,Ba,Tg,Kg,Bg,Tm2a,Bm,Vm,mag_strength]=ImuCalibration_Gesture_nplot(data,threshold);

% input data raw IMU data from mpu9250 
% data :time accelerometer  gyroscope   magnetometer 
%  cal_acc=Ta*Ka*(raw_acc+Ba)
%  cal_gyro=Tg*Kg*(raw_gyro+Bg)
%  cal_mag=Tm2a*(raw_mag+Bm)

%%% corr_data print
% text={'Ba' '' '' 'Ta' '' '' 'Ka'};
% name='corr_data';
% xlswrite('test.xlsx',text,name,'A1')
% xlswrite('test.xlsx',Ba,name,'A2')
% xlswrite('test.xlsx',Ta,name,'D2')
% xlswrite('test.xlsx',Ka,name,'G2')
% text={'Bg' '' '' 'Tg' '' '' 'Kg'};
% xlswrite('test.xlsx',text,name,'A5')
% xlswrite('test.xlsx',Bg,name,'A6')
% xlswrite('test.xlsx',Tg,name,'D6')
% xlswrite('test.xlsx',Kg,name,'G6')
% text={'Bm' '' '' 'Tm2a'};
% xlswrite('test.xlsx',text,name,'A9')
% xlswrite('test.xlsx',Bm,name,'A10')
% xlswrite('test.xlsx',Tm2a,name,'D10')

[a1 a2]=size(data);
data_corr=zeros(a1,a2);
data_corr(:,1)=data(:,1);
for i=1:a1   %校正數據
    data_corr(i,2:4)=(Ta*Ka*(data(i,2:4)'+Ba))';
    data_corr(i,5:7)=(Tg*Kg*(data(i,5:7)'+Bg))';
    data_corr(i,8:10)=(Tm2a*(data(i,8:10)'+Bm))';
end
X=data(:,8:10); 
Xnew=data_corr(:,8:10);
Xneww=sqrt(Xnew(:,1).^2+Xnew(:,2).^2+Xnew(:,3).^2);
Xw=sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2);
%%
close all
figure(1)
hold on
% 
scatter3(X(:,1),X(:,2),X(:,3),'red.')
scatter3(Xnew(:,1),Xnew(:,2),Xnew(:,3),'blue.')
scatter3(Bm(1)*-1,Bm(2)*-1,Bm(3)*-1,'ro')
scatter3(0,0,0,'blueo')
% 
title('橢球校正')
grid on
grid minor
axis equal tight
xlabel('X');
ylabel('Y');
zlabel('Z');
hold off
% 
figure(2)
hold on
plot(Xneww)
plot(Xw)
grid on
hold off
% %% data print
% 
% 
% [a1 a2]=size(data);
% data_corr=zeros(a1,a2);
% data_corr(:,1)=data(:,1);
% for i=1:a1
%     data_corr(i,2:4)=(Ta*Ka*(data(i,2:4)'+Ba))';
%     data_corr(i,5:7)=(Tg*Kg*(data(i,5:7)'+Bg))';
%     data_corr(i,8:10)=(Tm2a*(data(i,8:10)'+Bm))';
% end
% delemagn='data_corr';
% delemagdata=data_corr;
% text={'time' 'ax' 'ay' 'az' 'gx' 'gy' 'gz' 'mx' 'my' 'mz'};
% xlswrite('data_corr.xlsx',text,delemagn)
% xlswrite('data_corr.xlsx',delemagdata,delemagn,'A2')

%% old Corr_Mag

close all
% mo=load('C:\Users\henry chen\Desktop\LOG24.csv');
mo=load('proj\11 19 power_v4 box done test\57 54 48 36 24 12 6.csv');
% mo=load('proj\11 10 11 12 big c_v2 200v meeting\two big c compare 12 24 36 48 54_2.csv');
% mo=xlsread('proj\7 2 中油 iron dic and pull\點8 blank 2')

% mo1=load('proj\9 24 中油 iron pull big c\log0.csv');
% mo2=load('proj\9 24 中油 iron pull big c\log1.csv');
% mo3=load('proj\9 24 中油 iron pull big c\log2.csv');
% mo4=load('proj\9 24 中油 iron pull big c\log3.csv');
% mo5=load('proj\9 24 中油 iron pull big c\log4.csv');

% mo1=load('C:\Users\henry chen\Desktop\magtest\proj\9 24 中油 iron pull big c\log0.csv');
% mo2=load('C:\Users\henry chen\Desktop\magtest\proj\9 24 中油 iron pull big c\log1.csv');
% mo3=load('C:\Users\henry chen\Desktop\magtest\proj\9 24 中油 iron pull big c\log2.csv');
% mo4=load('C:\Users\henry chen\Desktop\magtest\proj\9 24 中油 iron pull big c\log3.csv');
% mo5=load('C:\Users\henry chen\Desktop\magtest\proj\9 24 中油 iron pull big c\log4.csv');

% mo=[mo1;mo2];
% mo=[mo3;mo4;mo5];
moo=mo;
m21=moo(:,11:13);
m22=moo(:,14:16);

% moo=load('C:\Users\henry chen\Desktop\magtest\proj\10 14 2mag big c\log4.csv');

[No,~]=size(moo); 
t0=[1 No];
smk=25; %移動平均係數
smk=1; %移動平均係數
testn='1'; %檔名
t=moo(:,1);

mo=moo(:,[1:4,8:13]);
[a1 a2]=size(mo);
mo_cor=zeros(a1,a2);
mo_cor(:,1)=mo(:,1);
for i=1:a1   %校正數據
    mo_cor(i,2:4)=(Ta*Ka*(mo(i,2:4)'+Ba))';
    mo_cor(i,5:7)=(Tg*Kg*(mo(i,5:7)'+Bg))';
    mo_cor(i,8:10)=(Tm2a*(mo(i,8:10)'+Bm))';
end

acc1 = mo_cor(:,2:4);
% acc2 = mo_cor(:,5:7);
m= mo_cor(:,8:10);

ax=acc1(:,1);ay=acc1(:,2);az=acc1(:,3);
% ax2=acc2(:,1);ay2=acc2(:,2);az2=acc2(:,3);
aGo=sqrt(mo_cor(:,2).^2+mo_cor(:,3).^2+mo_cor(:,4).^2);
% aG2o=sqrt(mo_cor(:,5).^2+mo_cor(:,6).^2+mo_cor(:,7).^2);
aG=sqrt(ax.^2+ay.^2+az.^2);
% aG2=sqrt(ax2.^2+ay2.^2+az2.^2);

GG=sqrt(mo_cor(:,5).^2+mo_cor(:,6).^2+mo_cor(:,7).^2)*0.01;

mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);

x=smooth(m(:,1),smk);
y=smooth(m(:,2),smk);
z=smooth(m(:,3),smk);
smedm=[x y z];
xyz=sqrt(x.^2+y.^2+z.^2);
%% differend pull in big data(but separated)
mooton=[];
mootoff=[];
% data=mo1;
% mooton=[8.8,11.8,31.6,35.7,55,58.1,63.9,67,73.8,87.4,109.6]*10^4; %mo1 %9 24 中油 iron pull
% mootoff=[10,13,33,36.6,56,59.15,64.8,67.8,74.6,88.5,110.8]*10^4;
% data=mo2;
% mooton=[2 4.4 8.6 11.3 14.3 17.2 19.9 22.4 25.4]*10^4; %mo2
% mootoff=[2.9 5.4 9.45 12.3 15.4 18.2 21 23.4 26.5]*10^4;
% data=mo3;
% mooton=[1,1.22,1.42,1.65,1.87,2.16,3.2,3.42,3.67,3.92,4.16,4.42,6.05,6.32,6.53,6.75]*10^5; %mo3
% mootoff=[1.09,1.3,1.52,1.73,1.97,2.26,3.3,3.5,3.76,4.01,4.26,4.53,6.13,6.41,6.63,6.85]*10^5;
% data=mo4;
% mooton=[5.22,5.44,6.36,6.58,6.77,7.06,7.5,8.13,10,10.57,10.9]*10^5; %mo4
% mootoff=[5.32,5.52,6.48,6.68,6.86,7.16,7.62,8.22,10.1,10.66,11]*10^5;
data=mo5;
mooton=[3.4 6.3 9.6 16 28 31.1 35.3 62 64.8 67.8]*10^4; %mo5
mootoff=[4.2 7.1 10.4 17.2 28.8 32 36.2 62.9 65.8 68.8]*10^4;

mooton=round(mooton);
mootoff=round(mootoff);

mot=[];
[N,k]=size(mooton);
for i=1:k
mot(i,1)=data(mooton(i),1);
mot(i,2)=data(mootoff(i),1);
end

pton=[];
ptoff=[];
for i=1:k
[~,pton(i)]=min(abs(t-mot(i,1)));
[~,ptoff(i)]=min(abs(t-mot(i,2)));
end

% ticks=[mot(:,1) mot(:,2)];
ticks=[pton ptoff];
ticks=reshape (ticks, 1, numel(ticks));
ticks=sort(ticks,'ascend');

% plot(t,m(:,1))
plot(m(:,1))
set(gca, 'xtick', ticks);
grid on
%% output data(excel)
filename='5';
for i=1:k %%
pulln=int2str(i);
% text={'time' 'a1x' 'a1y' 'a1z' 'a2x' 'a2y' 'a2z' 'gx' 'gy' 'gz' 'mx' 'my' 'mz' 'encoder' 'encodero1' 'encodero2'};
pulldata=[t(pton(i):ptoff(i)) mo_cor((pton(i):ptoff(i)),5:7) mo_cor((pton(i):ptoff(i)),2:4) mo_cor((pton(i):ptoff(i)),8:10)];
text={'time' 'ax' 'ay' 'az' 'gx' 'gy' 'gz' 'mx' 'my' 'mz'};
% xlswrite(filename,text,pulln)
xlswrite(filename,pulldata,pulln,'A2')
end
%%
for i=1:k %%
filename=[int2str(i),'.txt']
fid=fopen(['C:\Users\henry chen\Desktop\fly\Fly_IMUCalibration-Gesture-master\IMUCalibration-Gesture-master\data\',char(filename)],'w');%寫入檔案路徑
pulldata=[t(pton(i):ptoff(i)) mo_cor((pton(i):ptoff(i)),8:10) mo_cor((pton(i):ptoff(i)),2:4) mo_cor((pton(i):ptoff(i)),5:7)];
[r,c]=size(pulldata);            % 得到矩陣的行數和列數
 for w=1:r
  for v=1:c
  fprintf(fid,'%f\t',pulldata(w,v));
  end
  fprintf(fid,'\r\n');
 end
fclose(fid);
end
%% threshold filter waves for 2 set compare
close all
T1=[28000 197000];
T0=[33000 204000];
[~,k]=size(T0);

waveN=15;
magmin(waveN,k-1)=0;
magmax(waveN,k-1)=0;

useddata=z; % 用來判斷的數據

hold on
plot((aG(t0(1):t0(2))-0.5))
plot(useddata)
filtCutOff = 1; % 高通判斷數據使可判斷波峰谷
[b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
useddataL = filtfilt(b, a, useddata);
% useddataL = useddata; % 高速升降無須低通
plot(useddataL)
ff=0;
ff=7; %繼電器開關
for i=1:k %找出每段的波峰波谷
testdata=useddataL(T1(i):T0(i));
[~,kt]=size(testdata);
IndMin=find(diff(sign(diff(testdata)))>0)+ff;  %獲得局部最小值的位置(可微調位置
IndMax=find(diff(sign(diff(testdata)))<0)+ff;  %獲得局部最大值的位置

IndMin(waveN+1:size(IndMin))=[]; % 使同長度
IndMax(waveN+1:size(IndMax))=[];

% magmin(:,i)=IndMin+T1(i)-2; % 電容放電time fix 其他0
magmin(:,i)=IndMin+T1(i); % 彌補時間差
magmax(:,i)=IndMax+T1(i);

plot(magmin(:,i),useddata(magmin(:,i)),'r^') % 標示在原數據
plot(magmax(:,i),useddata(magmax(:,i)),'k*')
end
% plot(endcoder)

tag=[T0  T1]; % 取樣區間
% tag=[magmax(1,:) magmin(1,:)]; % 檢查位置
tag=reshape (tag, 1, numel(tag));
tag=sort(tag,'ascend');
set(gca, 'xtick', tag);
grid on
hold off
xlabel('Time')
ylabel('acc、LPF-mag、mag')


ew(k,:)=0;
for i=1:k %每個位置磁力平均值
    for j=1:waveN
        wx(j,i)=x(magmax(j,i))-x(magmin(j,i));
        wy(j,i)=y(magmax(j,i))-y(magmin(j,i));
        wz(j,i)=z(magmax(j,i))-z(magmin(j,i));
    end
    ew(i)=sqrt(mean(wx(:,i))^2+mean(wy(:,i))^2+mean(wz(:,i))^2);
    ews(i,:)=[mean(wx(:,i)) mean(wy(:,i)) mean(wz(:,i))];
end
ew;
%% 取得拉動平均速度
smf=300; %調整smooth參數判斷 
% aGsq=abs(aG-mode(aG(420000:500000))); %mo1~2
aGsq=abs(aG-mode(aG(800000:1000000))); %mo3`5
aGsq=smooth(aGsq(t0(1):t0(2)),smf);
stationary=aGsq>1.2; %設定閥值 使0為靜止 1為移動 
at=[];

close  all
hold on
plot(stationary)
plot(aGsq*0.1)
plot(mm4*0.01)
at(k,2)=0;
for i=1:k   
    at(i,1)=find(stationary(pton(i):ptoff(i)) > 0, 1, 'first')+pton(i);
    at(i,2)=find(stationary((pton(i)+ptoff(i))/2:ptoff(i)) == 0, 1, 'first')+(pton(i)+ptoff(i))/2;
end
at=round(at);

% atticks=[at(:,1),at(:,2),pton',ptoff'];
atticks=[at(:,1),at(:,2)];
atticks=reshape (atticks, 1, numel(atticks));
atticks=sort(atticks,'ascend');
set(gca, 'xtick', atticks);
grid on
hold off

D=9.5; %已知移動距離(m)
pv=zeros(k,1);
 for i=1:k %  距離/移動時間
     pv(i,1)=D/(t(at(i,2))-t(at(i,1)));
 end
 pv
 %% 2mag
 m2=moo(:,14:16);
for i=1:a1   %校正數據
    m2(i,1:3)=(Tm2a*(moo(i,14:16)'+Bm))';
end
mm42=sqrt(m2(:,1).^2+m2(:,2).^2+m2(:,3).^2);
aaxis=3;
% m2(:,aaxis)=m2(:,aaxis)+29.05;
figure(1)
close all
hold on
plot(m(:,aaxis))
plot(m2(:,aaxis))
hold off
figure(2)
hold on
plot(mm4)
plot(mm42)
hold off
figure(3)
hold on
plot(m)
plot(m2)
hold off
%%
filename=['2magrun','.txt']
fid=fopen(['data\',char(filename)],'w');%寫入檔案路徑
pulldata=[m m2];
[r,c]=size(pulldata);            % 得到矩陣的行數和列數
 for w=1:r
  for v=1:c
  fprintf(fid,'%f\t',pulldata(w,v));
  end
  fprintf(fid,'\r\n');
 end
fclose(fid);