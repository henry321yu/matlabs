clear; close all;


fitdata=xlsread('fit data\11 26 new fit data 2mag new place_th_2.8_m2濾y.xlsx');
Ba=fitdata(1:3,1);Ta=fitdata(1:3,2:4);Ka=fitdata(1:3,5:7);
Bg=fitdata(5:7,1);Tg=fitdata(5:7,2:4);Kg=fitdata(5:7,5:7);
Bm=fitdata(9:11,1);Tm2a=fitdata(9:11,2:4);
Bm_2=fitdata(13:15,1);Tm2a_2=fitdata(13:15,2:4);

% 時間  加速度(m/s/s)XYZ 角速度XYZ 磁力XYZ 放在data.txt 中
data=load('fit data\11 26 new fit code 2mag new place grass.csv');

acc=data(:,2:4);
acc=data(:,5:7);
gyr=data(:,8:10);
mag=data(:,11:13);
mag2=data(:,14:16);
Aa=sqrt(acc(:,1).^2+acc(:,2).^2+acc(:,3).^2);
Ag=sqrt(gyr(:,1).^2+gyr(:,2).^2+gyr(:,3).^2);
Am=sqrt(mag(:,1).^2+mag(:,2).^2+mag(:,3).^2);
Am2=sqrt(mag2(:,1).^2+mag2(:,2).^2+mag2(:,3).^2);


[a1 a2]=size(data);
data_corr=zeros(a1,a2);
data_corr(:,1)=data(:,1);
for i=1:a1   %校正數據
    data_corr(i,2:4)=(Ta*Ka*(data(i,2:4)'+Ba))';
    data_corr(i,8:10)=(Tg*Kg*(data(i,8:10)'+Bg))';
    data_corr(i,11:13)=(Tm2a*(data(i,11:13)'+Bm))';
    data_corr(i,14:16)=(Tm2a_2*(data(i,14:16)'+Bm_2))';
end
X=data(:,11:13); 
Xnew=data_corr(:,11:13);
Xneww=sqrt(Xnew(:,1).^2+Xnew(:,2).^2+Xnew(:,3).^2);
Xw=sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2);
X2=data(:,14:16); 
Xnew2=data_corr(:,14:16);
Xneww2=sqrt(Xnew2(:,1).^2+Xnew2(:,2).^2+Xnew2(:,3).^2);
Xw2=sqrt(X2(:,1).^2+X2(:,2).^2+X2(:,3).^2);

close all
figure(1)
hold on

scatter3(X(:,1),X(:,2),X(:,3),'red.')
scatter3(Xnew(:,1),Xnew(:,2),Xnew(:,3),'blue.')
scatter3(Bm(1)*-1,Bm(2)*-1,Bm(3)*-1,'ro')
scatter3(0,0,0,'blueo')

title('橢球校正X')
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
plot(Xw)
plot(Xneww)
grid on
hold off

figure(3)
hold on
% 
scatter3(X2(:,1),X2(:,2),X2(:,3),'red.')
scatter3(Xnew2(:,1),Xnew2(:,2),Xnew2(:,3),'blue.')
scatter3(Bm_2(1)*-1,Bm_2(2)*-1,Bm_2(3)*-1,'ro')
scatter3(0,0,0,'blueo')
% 
title('橢球校正X2')
grid on
grid minor
axis equal tight
xlabel('X');
ylabel('Y');
zlabel('Z');
hold off
% 
figure(4)
hold on
plot(Xw2)
plot(Xneww2)
grid on
hold off

figure(5)
hold on
plot(Xneww)
plot(Xneww2)
grid on
title('橢球校正X、X2 compare')
hold off

%% Corr_Mag
clear all;close all;

fitdata=xlsread('fit data\11 26 new fit data 2mag new place_th_2.8_m2濾y.xlsx');
Ba=fitdata(1:3,1);Ta=fitdata(1:3,2:4);Ka=fitdata(1:3,5:7);
Bg=fitdata(5:7,1);Tg=fitdata(5:7,2:4);Kg=fitdata(5:7,5:7);
Bm=fitdata(9:11,1);Tm2a=fitdata(9:11,2:4);
Bm_2=fitdata(13:15,1);Tm2a_2=fitdata(13:15,2:4);

close all
% mo1=load('proj\2020 proj\12 22 靜置 旋轉測陀螺儀\LOG64.csv');
% mo2=load('proj\2020 proj\12 22 靜置 旋轉測陀螺儀\LOG65.csv');
% mo3=load('proj\2020 proj\12 22 靜置 旋轉測陀螺儀\LOG66.csv');
% mo4=load('proj\2020 proj\12 22 靜置 旋轉測陀螺儀\LOG67.csv');
% mo5=load('proj\2020 proj\12 22 靜置 旋轉測陀螺儀\LOG68.csv');

% mo1=load('proj\5 25 26 羅森實測\5 25\log231.csv');
% mo2=load('proj\5 25 26 羅森實測\5 25\log232.csv');
% mo3=load('proj\5 25 26 羅森實測\5 25\log233.csv');
% mo4=load('proj\5 25 26 羅森實測\5 25\log234.csv');
% mo5=load('proj\5 25 26 羅森實測\5 25\log235.csv');
% mo6=load('proj\5 25 26 羅森實測\5 25\log236.csv');

% mo=load('C:\Users\henry chen\Desktop\temp\log343.csv');  %%%TEMP FOLDER%%%

% mo=load('proj\1 8 bigc_v2 test air and our pipe\12 24 36 48 54 54in iron.csv');
% mo=load('proj\1 7 中油鐵管 big c_v2\log77.csv');
% mo=load('proj\1 7 中油鐵管 big c_v2\log78_2.csv');
% mo=load('proj\1 7 中油鐵管 big c_v2\log82.csv');
% mo1=load('proj\1 7 中油鐵管 big c_v2\log83.csv');
% mo=load('proj\1 5 rosen itx504 d meeting room run\log73.csv');
% mo=load('proj\1 4 car run\office run2.csv');
% mo=load('proj\2020 proj\12 31 rosen itx504\504 in iron pipe 0.4 0.8 1.2 1.6  2 2.4 2.8 3.2.csv');
% mo=load('proj\5 4 MM\355 M 12 24 36 42 45 48 54 MM 12 24 36 42 45 48.csv');
% mo=load('proj\5 25 羅森實測\log168.csv');
% mo=load('proj\2020 proj\11 30 d to fly\log42.csv');
% mo=load('proj\2020 proj\11 26 fit and run\走廊boxpull.csv');
% mo=load('proj\2020 proj\11 19 power_v4 box done test\57 54 48 36 24 12 6.csv');
% mo=load('proj\2020 proj\11 10 11 12 big c_v2 200v meeting\two big c compare 12 24 36 48 54_2.csv');
% mo=xlsread('proj\7 2 中油 iron dic and pull\點8 blank 2')
% mo=xlsread('proj\6 25 giant wire peak for 5 26\LOG264.csv');
mo=xlsread('proj\6 28 giant wire d for 5 26\d for 5 26.csv');
% mo=xlsread('proj\6 28 giant wire d for 5 26\LOG265.csv');6 29 3 axis vector
% mo=xlsread('proj\6 29 3 axis vector\LOG270_.csv');
% mo=xlsread('proj\7 5 yz反轉嗎\遠 一群鐵 空氣.CSV');
% mo=xlsread('proj\7 7 hall test 10m\LOG286.CSV');
% mo=xlsread('proj\7 9 手推車 hall\LOG289.CSV');
% mo=xlsread('proj\7 12 time and sd speed\1號.CSV');
% mo=xlsread('proj\7 12 手推車hall M air iron\LOG7.CSV');

% moo=[mo1;mo2;mo3;mo4;mo5;mo6];

moo=mo;

[No,~]=size(moo); 
t0=[1 No];

% mo=moo(:,[1:4,8:16]);
mo=moo(:,[1:13]); %%16475
[a1 a3]=size(moo);
[a1 a2]=size(mo);
t=moo(:,1);  %time
f=moo(:,a3); %freqency
mo_cor=zeros(a1,a2);
mo_cor=mo;

mo_cor(:,1)=mo(:,1);
for i=1:a1   %校正數據
    mo_cor(i,2:4)=(Ta*Ka*(mo(i,2:4)'+Ba))';
    mo_cor(i,5:7)=(Tg*Kg*(mo(i,5:7)'+Bg))';
    mo_cor(i,8:10)=(Tm2a*(mo(i,8:10)'))';
    mo_cor(i,11:13)=(Tm2a_2*(mo(i,11:13)'+Bm_2))';
end

acc1 = mo_cor(:,2:4); %16475     已校正
m = mo_cor(:,8:10);
m2 = mo_cor(:,11:13);
gyr=mo_cor(:,5:7);

% acc1 = mo(:,2:4); %16475     無校正
% m = mo(:,8:10); %%%%
% m2 = mo(:,11:13); %%%%
% gyr=mo(:,5:7);



ax=acc1(:,1);ay=acc1(:,2);az=acc1(:,3);
% ax2=acc2(:,1);ay2=acc2(:,2);az2=acc2(:,3);
aGo=sqrt(mo_cor(:,2).^2+mo_cor(:,3).^2+mo_cor(:,4).^2);
% aG2o=sqrt(mo_cor(:,5).^2+mo_cor(:,6).^2+mo_cor(:,7).^2);
aG=sqrt(ax.^2+ay.^2+az.^2);
% aG2=sqrt(ax2.^2+ay2.^2+az2.^2);

GG=sqrt(gyr(:,1).^2+gyr(:,2).^2+gyr(:,3).^2)*0.01;

mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);
mm42=sqrt(m2(:,1).^2+m2(:,2).^2+m2(:,3).^2);

smk=1; %移動平均係數
x=smooth(m(:,1),smk);
y=smooth(m(:,2),smk);
z=smooth(m(:,3),smk);
smedm=[x y z];
xyz=sqrt(x.^2+y.^2+z.^2);
stdd=std(mo)';
stdd_cor=std(mo_cor)';
plot(mm4)
%% 7 12 time and sd speed
close all
clc

fc =20;fs=400; % Cut off frequency 7 12 3號
fc =50;fs=400; % Cut off frequency 7 12 3號
[j,p] = butter(2,fc/(fs/2),'low');

start=7650;n=500;
aa=start:start+n;

mo=xlsread('proj\7 12 time and sd speed\5 25.CSV');
t0=mo(aa,1);t0=t0-t0(1);
mo=xlsread('proj\7 12 time and sd speed\1號.CSV');
t1=mo(aa,1);t1=t1-t1(1);
mo=xlsread('proj\7 12 time and sd speed\3號.CSV');
t3=mo(aa,1);t3=t3-t3(1);
mo=xlsread('proj\7 12 time and sd speed\原 last.CSV');
t01=mo(aa,1);t01=t01-t01(1);
tt0=t0;tt1=t1;tt3=t3;tt01=t01;

[~,k]=size(aa);
for i=1:k-1
    tt0(i)=t0(i+1)-t0(i);
    tt1(i)=t1(i+1)-t1(i);
    tt3(i)=t3(i+1)-t3(i);
    tt01(i)=t01(i+1)-t01(i);
end
tt0(end-1:end)=[];tt1(end-1:end)=[];
tt3(end-1:end)=[];tt01(end-1:end)=[];
ltt0=filtfilt(j, p, tt0);
ltt1=filtfilt(j, p, tt1);
ltt3=filtfilt(j, p, tt3);
ltt01=filtfilt(j, p, tt01);

subplot(2,1,1)
hold on

plot(t1,aa,'.')
plot(t3,aa,'.')
plot(t01,aa,'.')
plot(t0,aa,'.')
hold off
legend('1','3','原 last','5 25')
xlabel('time')

subplot(2,1,2)
hold on
plot(tt1+0.03,'.')
plot(tt3+0.02,'.')
plot(tt01+0.01,'.')
plot(tt0,'.')
legend('1','3','原 last','5 25')
ylabel('每筆資料時間差(有加上偏移使分開)')

figure(2)
hold
plot(ltt1)  % LF
plot(ltt3)
plot(ltt01)
plot(ltt0)
legend('1','3','原 last','5 25')
ylabel('每筆資料時間差(LF)')
%% threshold filter (GGlow>0.005)(not use 因為還需要磁力對準)
aGs=mode(aG(1:20000));
movet=15000; %偵測到移動 下一次開始偵測時間
T0t=-5000; %偵測到移動 on(減)的時間
% T0t=5000; %偵測到移動 on(減)的時間
T1t=-500; %偵測到移動 off(加)的時間
T0=[];T1=[];
c=0;
fc = 1;fs=400; % Cut off frequency
[b,a] = butter(2,fc/(fs/2),'low');
GGlow=filtfilt(b, a, GG);
T0(end+1)=1;
for i=1:a1 %判斷移動後靜止的時間
    if GGlow(i)>0.005 && c==0
        T0(end+1)=i+T0t;
        T1(end+1)=i+T1t;
        c=1;
    end
    if GGlow(i)>0.005 && i>T0(end)+movet
        T0(end+1)=i+T0t;
        T1(end+1)=i+T1t;
    end
end
T1(end+1)=t0(2)-1;

% tticks=[T0 T1];
tticks=[T0];
tticks=reshape (tticks, 1, numel(tticks));
tticks=sort(tticks,'ascend');

close all
hold on
% plot(xyz)
% plot(xyz2)
plot(y)
plot(GG*100)
set(gca, 'xtick', tticks);
grid on
%%
close all
% hold on
% plot(aG*10.^2)
% plot(t,GG*10.^2)
% % plot(t,mm4)
% plot(t,x)
% plot(t,y)
% plot(t,z)
% plot(t,mo_cor(:,9)-Bm(1))
% plot(t,mo(:,9))
%% D
on1=[56.721 101.085 134.412 179.272 213.356 256.46 288.54 322.88 347.94 375.753 417.861 ... %1m
    448.692 474.503 496.308 534.15 563.978 598.061 631.9 665.482 696.806 738.912 780.767 ... %2m
    895.552 938.659 972.994 1011.096 1041.168 1070.243 1098.813 1133.4 1158.46 1187.032 1214.603]; %3m %giant wire d for 5 26
on=on1;
off=on+4;

[~,k]=size(on);
ed=[];
eds=[];
ton=[];
toff=[];
ed(k,:)=0;
eds(k,3)=0;
ton(k,:)=0;
toff(k,:)=0;

for i=1:k
[~,ton(i)]=min(abs(t-on(i)));
[~,toff(i)]=min(abs(t-off(i)));
end

% ton=[66499 179067 521486];
% toff=ton+1000;
% k=3

ticks=[ton toff];
ticks=reshape (ticks, 1, numel(ticks));
ticks=sort(ticks,'ascend');

for i=1:k % use time no filter
    ed(i)=sqrt((m(ton(i),1)-m(toff(i),1))^2+(m(ton(i),2)-m(toff(i),2))^2+(m(ton(i),3)-m(toff(i),3))^2);
    eds(i,:)=[m(ton(i),1)-m(toff(i),1) m(ton(i),2)-m(toff(i),2) m(ton(i),3)-m(toff(i),3)]';
end

ed
% eds

onoffacc1(k,3)=0; 
for i=1:k  %每個位置加速度
        onoffacc1(i,:)=mean(acc1(ton(i):toff(i),:));
end
onoffaccdata=onoffacc1;

close all
hold on
plot(m(:,1))
plot(m(:,2))
plot(m(:,3))

% plot(G)
% plot(aG)

tickson=[ton];
ticksoff=[toff];
tickson=reshape (tickson, 1, numel(tickson));
ticksoff=reshape (ticksoff, 1, numel(ticksoff));
tickson=sort(tickson,'ascend');
ticksoff=sort(ticksoff,'ascend');

for i=1:k
    plot(tickson(i),x(tickson(i)),'bo')
    plot(ticksoff(i),x(ticksoff(i)),'ro')
    plot(tickson(i),y(tickson(i)),'bo')
    plot(ticksoff(i),y(ticksoff(i)),'ro')
    plot(tickson(i),z(tickson(i)),'bo')
    plot(ticksoff(i),z(ticksoff(i)),'ro')
end
set(gca, 'xtick', ticks);
grid on
xlabel('Time')
ylabel('Gauss')
hold off

dd=[503.2000099 471.4318615 441.4612101 413.1924491 386.998708 363.0513049 342.2352992 325.3382855 312.5447968 304.7035444 302 ...
    201 205.0396303 216.5207842 234.6103152 257.5305807 284.610699 314.5870309 346.3018914 379.5853 414.0591745 449.8969326 ...
    414.736362 375.5582511 337.17058 299.2056149 261.8472837 224.9494388 189.528362 156.9745202 128.375426 107.8899903 100];  %giant wire d for 5 26

figure(2)
loglog(dd(1:11),ed(1:11))
figure(3)
loglog(dd(12:22),ed(12:22))
figure(4)
loglog(dd(23:33),ed(23:33))
figure(5)
loglog(dd,ed,'.')
close all
%%
% x=[]
% fun = @(x,xdata)x(1)*exp(x(2)*xdata);
% x0 = [700,5];
% x = lsqcurvefit(fun,x0,dd(12:22)',ed(12:22))




rng default % for reproducibility
xdata = linspace(0,3);
ydata = exp(-1.3*xdata) + 0.05*randn(size(xdata));
lb = [0,-2];
ub = [3/4,-1];
x0 = [1/2,-2];
x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub)
plot(xdata,ydata,'ko',xdata,fun(x,xdata),'b-')
legend('Data','Fitted exponential')
title('Data and Fitted Curve')

%% on off point
% start_n=[0.1 0.215 0.13 0.24 0.25]; %2m
% start_n=[0.103 0 0.131 0.054 0.197]; %3m
start_n=[0.071 0.188 0.131 0.12 0.21]; %4.5m
% start_n=0
Duration_t=0.2507;
round_N=floor((t(toff)-t(ton))/Duration_t)-1;
take=1
pulln=5

takeont=[];
takeon=[];
takeont(pulln,:)=0;
takeon(pulln,:)=0;

for i=1:round_N(take)
    takeont(i)=t(ton(take))+Duration_t*i;
end
takeont=takeont+start_n(take);
% takeont=takeont+start_n;

for i=1:round_N(take)
[~,takeon(i)]=min(abs(t-takeont(i)));
end

% close all
hold on
% plot(m(ton(take):toff(take),1))
% % plot(m(ton(take):toff(take),2))
% % plot(m(ton(take):toff(take),3))
% set(gca, 'xtick', takeon-ton(take));

yzironbia=3

plot(t(ton(take):toff(take)),m(ton(take):toff(take),1))
plot(t(ton(take):toff(take)),m(ton(take):toff(take),2))
plot(t(ton(take):toff(take)),m(ton(take):toff(take),3))
for i=1:round_N(take)
plot(t(takeon(i)),m(takeon(i),1),'r.')
plot(t(takeon(i)-yzironbia),m(takeon(i)-yzironbia,2),'b.')
plot(t(takeon(i)-yzironbia),m(takeon(i)-yzironbia,3),'g.')
end

% plot(ton(take):toff(take),m(ton(take):toff(take),1)) % for search bia
% for i=1:round_N(take)
% plot(takeon(i),m(takeon(i),1),'r.')
% end

set(gca, 'xtick', takeont);
% set(gca,'xTickLabel',num2str(get(gca,'xTick')','%.2f'))
grid on
xlabel('time(s)')
ylabel('mag(uT)')
hold off
%% threshold filter waves for 2 set compare
close all
T1=ton;
T0=toff;
[k,~]=size(T0);
waveN=7 %計算的波數
magmin=[];
magmax=[];
magmin(waveN,k-1)=0;
magmax(waveN,k-1)=0;

useddata=x; % 用來判斷的數據

hold on
plot((aG(t0(1):t0(2))-0.5))
% plot(useddata)
filtCutOff = 1; % 高通判斷數據使可判斷波峰谷
[b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
useddataL = filtfilt(b, a, useddata);
% plot(useddataL)
ff=-2;
for i=1:k %找出每段的波峰波谷
testdata=useddataL(T1(i):T0(i));
[~,kt]=size(testdata);
IndMin=find(diff(sign(diff(testdata)))>0)+ff;  %獲得局部最小值的位置(可微調位置
IndMax=find(diff(sign(diff(testdata)))<0)+ff;  %獲得局部最大值的位置

IndMin(waveN+1:size(IndMin))=[]; % 使同長度
IndMax(waveN+1:size(IndMax))=[];

magmin(:,i)=IndMin+T1(i); % 彌補時間差
magmax(:,i)=IndMax+T1(i);

plot(magmin(:,i),useddata(magmin(:,i)),'r^') % 標示在原數據
plot(magmax(:,i),useddata(magmax(:,i)),'k*')
end
% plot(endcoder)

plot(m)


% tag=[T0  T1]; % 取樣區間
tag=[magmax magmin]; % 檢查位置
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
%         wy(j,i)=y(magmax(j,i))-y(magmin(j,i));
%         wz(j,i)=z(magmax(j,i))-z(magmin(j,i));
        wy(j,i)=y(magmax(j,i)-yzironbia)-y(magmin(j,i)-yzironbia);
        wz(j,i)=z(magmax(j,i)-yzironbia)-z(magmin(j,i)-yzironbia);
    end
    ew(i)=sqrt(mean(wx(:,i))^2+mean(wy(:,i))^2+mean(wz(:,i))^2);
    ews(i,:)=[mean(wx(:,i)) mean(wy(:,i)) mean(wz(:,i))];
end
ews(1:14,:)=ews(1:14,:)*-1; %%%因參考值波峰波谷反轉故*-1
figure(2)
plot(ews)

onoffacc1(k,3)=0; 
for i=1:k  %每個位置加速度
        onoffacc1(i,:)=mean(acc1(ton(i):toff(i),:));
end
%%
smf=50; %調整smooth參數判斷 
% aGsq=abs(aG-mode(aG(15000:35000))); %2m
aGsq=abs(aG-mode(aG(1:200000))); %3m、4m
stationary=aGsq>0.0014; %設定閥值 使0為靜止 1為移動 
at=[];
pton=ton;
ptoff=toff;

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

D=11.38; %已知移動距離(m)
pv=zeros(k,1);
 for i=1:k %  距離/移動時間
     pv(i,1)=D/(t(at(i,2))-t(at(i,1)));
 end
 pv
 
 
 
 
 
 
 
 
 
 
 %% hall %%下多似計步器推算
 addpath('hall')
H=moo; %因以加速度濾波及編碼器計算 須修正至量測區間
% H=moo(35000:42000,:); %因以加速度濾波及編碼器計算 須修正至量測區間
H=moo(80000:87000,:); %因以加速度濾波及編碼器計算 須修正至量測區間
[Nf,~]=size(H);
h1=H(:,14);
h2=H(:,15);

% fc =10;fs=400; % Cut off frequency
% fc =10;fs=350; % Cut off frequency 7 12
fc =10;fs=340; % Cut off frequency 7 12 3號
[b,a] = butter(2,fc/(fs/2),'low');
% h1=filtfilt(b, a, h1);

map=[-10 10]; %map係數
mh1=(h1-min(h1))*(map(2)-map(1))/(max(h1)-min(h1))+map(1);
mh2=(h2-min(h2))*(map(2)-map(1))/(max(h2)-min(h2))+map(1);

close all
hold on
% plot(mh1)
plot(mh2)
% plot(mm4/100)
% plot(aG)
t=H(:,1);
tt=t-t(1);

inputEncoder1 = [tt mh1];
inputEncoder2 = [tt mh2];
Encoder_threshold=11; %%%% 參數須設定
close all
[Encoder_X1,turning_velocity1,turning_degree1] = Encoder_1(inputEncoder1,Encoder_threshold);

close all
[Encoder_X2,turning_velocity2,turning_degree2] = Encoder_1(inputEncoder2,Encoder_threshold);
Encoder_X1(Nf);
Encoder_X2(Nf)
endcoder=Encoder_X1;
if endcoder(Nf)<Encoder_X2(Nf)
    endcoder=Encoder_X2;
end
% endcoder=smooth(Encoder_X2,200);  %%%%%%%%%%%%%
%%
close  all

% figure(1)
% hold on
% % plot(tt,Encoder_X1);
% plot(tt,Encoder_X2);
% legend('hall1','hall2','Location','southeast')
% title('Encoder')
% hold off

figure(2)
hold on
plot(endcoder,smooth(H(:,11),25)); %m的x軸
plot(endcoder,smooth(H(:,12),25)); %m的x軸
plot(endcoder,smooth(H(:,13),25)); %m的x軸
xlabel('postion(m) (hall2)')
ylabel('magnetic field strength(guess) (x)')
legend('mag1','mag6','Location','southeast')
title('magnetic field strength in each position')
hold off

% stopt=[1 30];  %% 設定靜止時間(資料數)
% dt=1/fs;
% [M_velo2, dis_from_Mag2]=getVtromMag(m2*-1, m, dt, No ,stopt);
% M_dis2=cumtrapz(tt,M_velo2);
% M_dis2(No);   %  getVtromMag  積分出來的距離

% figure(3) %hall距離 vs 磁場推算距離
% hold on
% plot(tt,Encoder_X2);
% plot(tt,M_dis2);
% xlabel('time(s)')
% ylabel('postion(m)  (hall2)')
% legend('hall2 sensor','getVfromMag','Location','southeast')
% title('Encoder position and getVfromMag position')
% hold off

% figure(4) %hall距離 vs 磁場推算距離
% hold on
% plot(tt,turning_velocity1);
% plot(tt,M_velo2);
% xlabel('time(s)')
% ylabel('postion(m)  (hall2)')
% legend('hall2 sensor','getVfromMag','Location','southeast')
% title('Encoder velocity and getVfromMag velocity')
% hold off