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

% mo=load('proj\1 8 bigc_v2 test air and our pipe\12 24 36 48 54 54in iron.csv');
% mo=load('proj\1 7 中油鐵管 big c_v2\log77.csv');
% mo=load('proj\1 7 中油鐵管 big c_v2\log78_2.csv');
% mo=load('proj\1 7 中油鐵管 big c_v2\log82.csv');
% mo1=load('proj\1 7 中油鐵管 big c_v2\log83.csv');
% mo=load('proj\1 5 rosen itx504 d meeting room run\log73.csv');
% mo=load('proj\1 4 car run\office run2.csv');
% mo=load('proj\2020 proj\12 31 rosen itx504\504 in iron pipe 0.4 0.8 1.2 1.6  2 2.4 2.8 3.2.csv');
% mo=load('C:\Users\henry chen\Desktop\temp\log168.csv');
% mo=load('proj\5 4 MM\355 M 12 24 36 42 45 48 54 MM 12 24 36 42 45 48.csv');
% mo=load('proj\5 25 羅森實測\log168.csv');
% mo=load('proj\2020 proj\11 30 d to fly\log42.csv');
% mo=load('proj\2020 proj\11 26 fit and run\走廊boxpull.csv');
% mo=load('proj\2020 proj\11 19 power_v4 box done test\57 54 48 36 24 12 6.csv');
% mo=load('proj\2020 proj\11 10 11 12 big c_v2 200v meeting\two big c compare 12 24 36 48 54_2.csv');
% mo=xlsread('proj\7 2 中油 iron dic and pull\點8 blank 2')
% mo=xlsread('proj\6 25 giant wire peak for 5 26\LOG264.csv');
mo=xlsread('proj\6 28 giant wire d for 5 26\LOG267.csv');


% for i=1:100
%     mo(end,:)=[];
% end
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

% mo_cor(:,1)=mo(:,1);
% for i=1:a1   %校正數據
%     mo_cor(i,2:4)=(Ta*Ka*(mo(i,2:4)'+Ba))';
%     mo_cor(i,5:7)=(Tg*Kg*(mo(i,5:7)'+Bg))';
%     mo_cor(i,8:10)=(Tm2a*(mo(i,8:10)'+Bm))';
%     mo_cor(i,11:13)=(Tm2a_2*(mo(i,11:13)'+Bm_2))';
% end

% acc1 = mo_cor(:,2:4);
% acc2 = mo_cor(:,5:7);
% m = mo_cor(:,8:10);
% m2 = mo_cor(:,11:13);
% gyr=mo_cor(:,5:7);

acc1 = mo(:,2:4); %16475
m = mo(:,8:10); %%%%
m2 = mo(:,11:13); %%%%
gyr=mo(:,5:7);

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

plot(mm42)
%%
% 確認軸向
% aaa=acc1(:,1);aaa2=moo(:,2)*9.8;
% aaa=acc1(:,2);aaa2=moo(:,3)*9.8;
% aaa=acc1(:,3);aaa2=moo(:,4)*9.8;
% aaa=gyr(:,1);aaa2=moo(:,8)/131;
% aaa=gyr(:,2);aaa2=moo(:,9)/131;
% aaa=gyr(:,3);aaa2=moo(:,10)/131;
% aaa=m(:,1);aaa2=moo(:,11);
% aaa=m(:,2);aaa2=moo(:,12);
% aaa=m(:,3);aaa2=moo(:,13);
% close all
% hold on
% plot(aaa)
% plot(aaa2)
% hold off

% 17632   %比較單個波形 %158
% a=mm4(ans:ans+80,:)
% 20443
% b=mm4(ans:ans+80,:)
% 
% plot(a)
% hold
% plot(b)

%%
close all
hold on
plot(mm42)
% plot(aG*10.^2)
plot(GG*10.^2)
%%
% biad=10680;biam=0.8; %78_2   %設定起始點、偏移量
biad=14272;biam=0.8; %82、83
% biad=2300;biam=1.1;
% biad=15470;biam=0.9;

biay=biad-floor(biad/1200)*1200;
fixy=[];
round_fix=floor(No/1200); %去m2(y)雜訊
fixy(round_fix,1)=0;
nm2=mo_cor(:,11:13);
for i=1:round_fix
    fixy(i)=biay+1200*(i-1);
end
for i=1:round_fix
    nm2(fixy(i),2)=mo_cor(fixy(i),12)+biam;
end
for i=1:round_fix
    nm2(fixy(i)+1,2)=mo_cor(fixy(i)+1,12)+biam;
end
close all
hold on
plot(mo_cor(:,12))
plot(nm2(:,2))
hold off
m2=nm2;
% m2=mo_cor(:,11:13);
%%
% mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);
% mm42=sqrt(m2(:,1).^2+m2(:,2).^2+m2(:,3).^2);
% mm4=sqrt(m(:,1).^2+m(:,3).^2);  %% 不看Y軸
% mm42=sqrt(m2(:,1).^2+m2(:,3).^2);
mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2); %手動校正
mm42=sqrt((m2(:,1)-2.66+0.98).^2+(m2(:,2)+11.575+0.942).^2+(m2(:,3)+3.08).^2);

close all

figure(1)
hold on
plot(t,m(:,1))
% plot(t,Tm2(:,1))
plot(t,m2(:,1)-2.66+0.98)

% plot(m(:,1))
% plot(m2(:,1))
hold off

figure(2)
hold on
plot(t,m(:,2))
% plot(t,m2(:,2))
plot(t,m2(:,2)+11.575+0.942)

% plot(m(:,2))
% plot(m2(:,2))
hold off

figure(3)
hold on
plot(t,m(:,3))
% % plot(t,m2(:,3))
plot(t,m2(:,3)+3.08)

% plot(m(:,3))
% plot(m2(:,3))
hold off

figure(4)
hold on
% plot(t,mm4)
% plot(t,mm42)
plot(mm4)
plot(mm42)
plot(aG)
% plot(t,GG*1000)
% plot(t,aaa2*10000)
hold off

% close all  %陀螺儀角度
% filtCutOff = 0.0000000035;
% filtCutOff = 0.00000000203694;
% filtCutOff = 0.0000000020369403;
% filtCutOff = 0.0000000020369403039;
% filtCutOff = 0.0000000020369403042;
% filtCutOff = 0.00000000203694030424; %辦公室 90 180....
% order=1;



% filtCutOff = 0.00000000203694030424;
% order=1;
% 
% % filtCutOff = 0.000000000012215;
% % order=2;
% [b, a] = butter(order, (2*filtCutOff)/(1/350), 'high');
% aaa = filtfilt(b, a, gyr(:,3));
% % aaa=gyr(:,3);
% aaa2=cumtrapz(t,aaa);
% % st=1;en=90000;
% % aaa2(en)*180/pi-aaa2(st)*180/pi
% aaa3=mod(aaa2*180/pi,360); %求餘數
% hold on
% % plot(aaa2(37000:45000)*180/pi)
% % plot(aaa2*180/pi)
% plot(aaa3)
% % plot(aG*100)
% hold off
%%
% plot(m)
% plot(m2)
% plot(acc1(:,1))
% hold off
% a
% close all
% figure(6)
% hold on
% plot(t,m)
% plot(t,m2)
% % plot(m)
% % plot(acc1(:,1))
% hold off
% 
% close all
% figure(7)
% hold on
% plot(m(87000:100000,1)+40)  % 1 8 bigc_v2 test
% plot(m(87000:100000,2))
% plot(m(87000:100000,3))
% 
% plot(m(55000:68000,1)-80)
% plot(m(55000:68000,2))
% plot(m(55000:68000,3))
% hold off
% 
% % delemagn=int2str(i);
% delemagdata=moo(1:350001,:);
% xlswrite('proj\1 7 中油鐵管 big c_v2\log78_2.xlsx',delemagdata)
% % writematrix(delemagdata,'proj\1 7 中油鐵管 big c_v2\log78_3.csv')
%% 1 7 2m D
mins=[0 1 1 1 2 2 2 3 3 4 4 4 4 5 5 5 6 6 7 8 9 9 10 11 11 12 13 13 14 14 15 15 15 16 16 16 17 17 17 18 18];
secs=[40 8 28 52 14 33 55 15 40 0 15 35 55 12 33 51 24 44 53 26 0 30 12 0 35 33 23 54 16 45 11 30 53 18 35 53 16 33 56 23 53];

[~,k]=size(mins);
on=[];
for i=1:k
on(i)=mins(i)*60+secs(i);
end

% on=[40 68 88 112 134 153 175 195 220 240 255 275 295 312 333 351 384 404 473 505 540 570 612 660 695 753 803 834 856 885 911 930 953 978 995 1013 1036 1053 1075 1103 1133];
on=on+221.4360+2;
off=on+2;
%% 1 7 pull
% mins=[5 6 8 8 10]; %2m
% secs=[32 53 0 48 35];
mins=[27 29 30 31 32];%3m
secs=[43 5 5 13 11];
% mins=[46 47 48 50 51]; %4.5m
% secs=[19 23 46 2 0];

[~,k]=size(mins);
on=[];
for i=1:k
on(i)=mins(i)*60+secs(i);
end

% on=on+0;  %2m
on=on+221.436; %3m、4m
off=on+15;
%%
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
ticks=[ton toff];
ticks=reshape (ticks, 1, numel(ticks));
ticks=sort(ticks,'ascend');

for i=1:k % use time no filter
    ed(i)=sqrt((m(ton(i),1)-m(toff(i),1))^2+(m(ton(i),2)-m(toff(i),2))^2+(m(ton(i),3)-m(toff(i),3))^2);
    eds(i,:)=[m(ton(i),1)-m(toff(i),1) m(ton(i),2)-m(toff(i),2) m(ton(i),3)-m(toff(i),3)]';
end

% ed
% eds

onoffacc1(k,3)=0; 
for i=1:k  %每個位置加速度
        onoffacc1(i,:)=mean(acc1(ton(i):toff(i),:));
end
onoffaccdata=onoffacc1;

close all
hold on
plot(m(:,1))
% plot(m(:,2))
% plot(m(:,3))

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
end
% set(gca, 'xtick', ticks);
% grid on
xlabel('Time')
ylabel('Gauss')
hold off
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
%% start_n to fly
for i=1:5

take=i;
    
takeont=[];
takeon=[];
start_n_to_fly(pulln,:)=0;
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

start_n_to_fly(take)=takeon(1)-ton(take)
% start_n_to_fly(take)=abs(takeon(1)-ton(take)-87)
end
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
%% 1 7 pull data (excel)
for i=1:k
% delemagn=(100+(i-1)*50)*0.1; %11 30 1m-0.5m-4m
delemagn=int2str(i);
delemagdata=[t(ton(i):toff(i),:) acc1(ton(i):toff(i),:) gyr(ton(i):toff(i),:) m(ton(i):toff(i),:) m2(ton(i):toff(i),:)];
text={'time' 'a1x' 'a1y' 'a1z' 'gx' 'gy' 'gz' 'm1x' 'm1y' 'm1z' 'm2x' 'm2y' 'm2z'};
xlswrite('proj\1 7 中油鐵管 big c_v2\1 7 output 4.5m.xlsx',text,delemagn)
xlswrite('proj\1 7 中油鐵管 big c_v2\1 7 output 4.5m.xlsx',delemagdata,delemagn,'A2')
end
%% 1 7 2m D data (excel)
% delemagn=(100+(i-1)*50)*0.1; %11 30 1m-0.5m-4m
% delemagn=int2str(delemagn);
delemagn='定點';
delemagdata=[(1:k)' onoffacc1 ews];
text={'distance' 'ax' 'ay' 'az' 'm1x' 'm1y' 'm1z'};
xlswrite('proj\1 7 中油鐵管 big c_v2\1 7 d w yzbia.xlsx',text,delemagn)
xlswrite('proj\1 7 中油鐵管 big c_v2\1 7 d w yzbia.xlsx',delemagdata,delemagn,'A2')
%% data (txt)
for i=1:k %%
filename=[int2str(i),'.txt']
fid=fopen(['proj\1 7 中油鐵管 big c_v2\txt\',char(filename)],'w');%寫入檔案路徑
pulldata=[t(ton(i):toff(i),:) acc1(ton(i):toff(i),:) gyr(ton(i):toff(i),:) m(ton(i):toff(i),:) m2(ton(i):toff(i),:)];
[r,c]=size(pulldata);            % 得到矩陣的行數和列數
 for w=1:r
  for v=1:c
  fprintf(fid,'%f\t',pulldata(w,v));
  end
  fprintf(fid,'\r\n');
 end
fclose(fid);
end
%% 取得拉動平均速度
smf=50; %調整smooth參數判斷 
% aGsq=abs(aG-mode(aG(15000:35000))); %2m
aGsq=abs(aG-mode(aG(50000:90000))); %3m、4m
aGsq=smooth(aGsq(t0(1):t0(2)),smf);
stationary=aGsq>1.2; %設定閥值 使0為靜止 1為移動 
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
%  %% 6 25 peak for 5 25.26
% delemagn='peak';
% delemagdata=[t m];
% % delemagdata=delemagdata(643799:643901,:);
% delemagdata=delemagdata(650131:650298,:);
% text={'time' 'm1x' 'm1y' 'm1z'};
% xlswrite('proj\6 25 giant wire peak for 5 26\peak.xlsx',text,delemagn,'E1')
% xlswrite('proj\6 25 giant wire peak for 5 26\peak.xlsx',delemagdata,delemagn,'E2')