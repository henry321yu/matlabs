%% 橢球校正
clear all

XL=load('fit data\4 15 proj new2b with hall without M.txt');
TWuT=45; % (高斯G)需設定 % 台灣磁場強度背景值45uT (0.45高斯) 
equals = ''; % no constraints by default


X=XL(:,[11 12 13]); % mag1
[center, radii, evecs, v, chi2 ,corr_center ,corr_radii] = ellipsoid_fit_new( X, equals );
Xnew=(X-corr_center)./corr_radii*TWuT;

corr_center_X=center(1,1); % 取得 
corr_center_Y=center(2,1);
corr_center_Z=center(3,1);
corr_radii_X=radii(1,1);
corr_radii_Y=radii(2,1);
corr_radii_Z=radii(3,1);


% X2=XL(:,[14 15 16]); % mag2
% [center, radii, evecs, v, chi2 ,corr_center ,corr_radii] = ellipsoid_fit_new( X2, equals );
% Xnew2=(X2-corr_center)./corr_radii*TWuT;
% 
% corr_center_X2=center(1,1);
% corr_center_Y2=center(2,1);
% corr_center_Z2=center(3,1);
% corr_radii_X2=radii(1,1);
% corr_radii_Y2=radii(2,1);
% corr_radii_Z2=radii(3,1);

%%%%%%%%%%% other 4mag
% 
% X3=XL(:,[17 18 19]); % mag3
% [center, radii, evecs, v, chi2 ,corr_center ,corr_radii] = ellipsoid_fit_new( X3, equals );
% Xnew3=(X3-corr_center)./corr_radii*TWuT;
% 
% corr_center_X3=center(1,1);
% corr_center_Y3=center(2,1);
% corr_center_Z3=center(3,1);
% corr_radii_X3=radii(1,1);
% corr_radii_Y3=radii(2,1);
% corr_radii_Z3=radii(3,1);
% 
% 
% X4=XL(:,[20 21 22]); % mag4
% [center, radii, evecs, v, chi2 ,corr_center ,corr_radii] = ellipsoid_fit_new( X4, equals );
% Xnew4=(X4-corr_center)./corr_radii*TWuT;
% 
% corr_center_X4=center(1,1);
% corr_center_Y4=center(2,1);
% corr_center_Z4=center(3,1);
% corr_radii_X4=radii(1,1);
% corr_radii_Y4=radii(2,1);
% corr_radii_Z4=radii(3,1);
% 
% 
% X5=XL(:,[23 24 25]); % mag5
% [center, radii, evecs, v, chi2 ,corr_center ,corr_radii] = ellipsoid_fit_new( X5, equals );
% Xnew5=(X5-corr_center)./corr_radii*TWuT;
% 
% corr_center_X5=center(1,1);
% corr_center_Y5=center(2,1);
% corr_center_Z5=center(3,1);
% corr_radii_X5=radii(1,1);
% corr_radii_Y5=radii(2,1);
% corr_radii_Z5=radii(3,1);
% 
% 
% X6=XL(:,[26 27 28]); % mag6
X6=XL(:,[11 12 13]); % mag6
[center, radii, evecs, v, chi2 ,corr_center ,corr_radii] = ellipsoid_fit_new( X6, equals );
Xnew6=(X6-corr_center)./corr_radii*TWuT;

corr_center_X6=center(1,1);
corr_center_Y6=center(2,1);
corr_center_Z6=center(3,1);
corr_radii_X6=radii(1,1);
corr_radii_Y6=radii(2,1);
corr_radii_Z6=radii(3,1);
%%% 校正圖
Xneww=sqrt(Xnew(:,1).^2+Xnew(:,2).^2+Xnew(:,3).^2);
% Xneww2=sqrt(Xnew2(:,1).^2+Xnew2(:,2).^2+Xnew2(:,3).^2);
% Xneww3=sqrt(Xnew3(:,1).^2+Xnew3(:,2).^2+Xnew3(:,3).^2);
% Xneww4=sqrt(Xnew4(:,1).^2+Xnew4(:,2).^2+Xnew4(:,3).^2);
% Xneww5=sqrt(Xnew5(:,1).^2+Xnew5(:,2).^2+Xnew5(:,3).^2);
% Xneww6=sqrt(Xnew6(:,1).^2+Xnew6(:,2).^2+Xnew6(:,3).^2);


figure(1)
hold on
% 
scatter3(corr_center_X,corr_center_Y,corr_center_Z,'ro')
% scatter3(corr_center_X2,corr_center_Y2,corr_center_Z2,'bo')
% scatter3(corr_center_X3,corr_center_Y3,corr_center_Z3,'go')
% scatter3(corr_center_X4,corr_center_Y4,corr_center_Z4,'yo')
% scatter3(corr_center_X5,corr_center_Y5,corr_center_Z5,'mo')
% scatter3(corr_center_X6,corr_center_Y6,corr_center_Z6,'co')
scatter3(X(:,1),X(:,2),X(:,3),'red.')
% scatter3(X2(:,1),X2(:,2),X2(:,3),'blue.')
% scatter3(X3(:,1),X3(:,2),X3(:,3),'green.')
% scatter3(X4(:,1),X4(:,2),X4(:,3),'yellow.')
% scatter3(X5(:,1),X5(:,2),X5(:,3),'m.')
% scatter3(X6(:,1),X6(:,2),X6(:,3),'c.')
scatter3(Xnew(:,1),Xnew(:,2),Xnew(:,3),'red.')
% scatter3(Xnew2(:,1),Xnew2(:,2),Xnew2(:,3),'blue.')
% scatter3(Xnew3(:,1),Xnew3(:,2),Xnew3(:,3),'green.')
% scatter3(Xnew4(:,1),Xnew4(:,2),Xnew4(:,3),'yellow.')
% scatter3(Xnew5(:,1),Xnew5(:,2),Xnew5(:,3),'m.')
% scatter3(Xnew6(:,1),Xnew6(:,2),Xnew6(:,3),'c.')
scatter3(0,0,0,'blacko')
% 
title('橢球校正')
grid on
grid minor
axis equal tight
xlabel('X');
ylabel('Y');
zlabel('Z');
% legend('1','2','3','4','5','6');
hold off
% 
figure(2)
hold on
plot(Xneww)
% plot(Xneww2)
% plot(Xneww3)
% plot(Xneww4)
% plot(Xneww5)
% plot(Xneww6)
% legend('1','2','3','4','5','6');
hold off

% for i=1:3:8257 %animetion
%     hold on
% scatter3(X(i,1),X(i,2),X(i,3),'red.')
% pause(0.001);
% end
%% read
close all
% mo=load('C:\Users\henry chen\Desktop\LOG20.txt');
mo=load('proj\8 4 end huge  wire d\air.txt');

[No,~]=size(mo); 
t0=[1 No];
smk=25; %移動平均係數
testn='1'; %檔名
t=mo(:,1);

% t0=[63000 204000];%定點1 %編碼器用
% t0=[30000 162000];%定點2
% t0=[24000 228000];%定點3
% t0=[45000 212000];%定點4

m=mo(:,[11 12 13]);%% 一號磁力儀
m(:,1) = (m(:,1) - corr_center_X )./ corr_radii_X *TWuT;
m(:,2) = (m(:,2) - corr_center_Y )./ corr_radii_Y *TWuT;
m(:,3) = (m(:,3) - corr_center_Z )./ corr_radii_Z *TWuT;
m=m*0.01; %to Guess
% 
% m2=mo(:,[14 15 16]);%% 二號磁力儀
% m2(:,1) = (m2(:,1) - corr_center_X2 )./ corr_radii_X2 *TWuT;
% m2(:,2) = (m2(:,2) - corr_center_Y2 )./ corr_radii_Y2 *TWuT;
% m2(:,3) = (m2(:,3) - corr_center_Z2 )./ corr_radii_Z2 *TWuT;
% m2=m2*0.01;
% 
% m3=mo(:,[17 18 19]);%% 三號磁力儀
% m3(:,1) = (m3(:,1) - corr_center_X3 )./ corr_radii_X3 *TWuT;
% m3(:,2) = (m3(:,2) - corr_center_Y3 )./ corr_radii_Y3 *TWuT;
% m3(:,3) = (m3(:,3) - corr_center_Z3 )./ corr_radii_Z3 *TWuT;
% m3=m3*0.01;
% 
% m4=mo(:,[20 21 22]);%% 四號磁力儀
% m4(:,1) = (m4(:,1) - corr_center_X4 )./ corr_radii_X4 *TWuT;
% m4(:,2) = (m4(:,2) - corr_center_Y4 )./ corr_radii_Y4 *TWuT;
% m4(:,3) = (m4(:,3) - corr_center_Z4 )./ corr_radii_Z4 *TWuT;
% m4=m4*0.01;
% 
% m5=mo(:,[23 24 25]);%% 五號磁力儀
% m5(:,1) = (m5(:,1) - corr_center_X5 )./ corr_radii_X5 *TWuT;
% m5(:,2) = (m5(:,2) - corr_center_Y5 )./ corr_radii_Y5 *TWuT;
% m5(:,3) = (m5(:,3) - corr_center_Z5 )./ corr_radii_Z5 *TWuT;
% m5=m5*0.01;
% 
% m6=mo(:,[26 27 28]);%% 六號磁力儀
% m6=mo(:,[14 15 16]);%% 六號磁力儀
% m6(:,1) = (m6(:,1) - corr_center_X6 )./ corr_radii_X6 *TWuT;
% m6(:,2) = (m6(:,2) - corr_center_Y6 )./ corr_radii_Y6 *TWuT;
% m6(:,3) = (m6(:,3) - corr_center_Z6 )./ corr_radii_Z6 *TWuT;
% m6=m6*0.01;

%%% acc corr fit
corr_center_aX=0.0171511511529850;
corr_center_aY=-0.00318848470804418;
corr_center_aZ=-0.0142065034430159;
corr_radii_aX=1.03431459019208;
corr_radii_aY=1.01712801230964;
corr_radii_aZ=1.01255711628185;
corr_center_aX2=0.00932117486074120;
corr_center_aY2=-0.00399760607468670;
corr_center_aZ2=0.00609017619411123;
corr_radii_aX2=1.01476072282633;
corr_radii_aY2=1.00927669187293;
corr_radii_aZ2=1.00009425512849;

acc1(:,1) = (mo(:,2) - corr_center_aX )./ corr_radii_aX;
acc1(:,2) = (mo(:,3) - corr_center_aY )./ corr_radii_aY;
acc1(:,3) = (mo(:,4) - corr_center_aZ )./ corr_radii_aZ;
acc2(:,1) = (mo(:,5) - corr_center_aX2 )./ corr_radii_aX2;
acc2(:,2) = (mo(:,6) - corr_center_aY2 )./ corr_radii_aY2;
acc2(:,3) = (mo(:,7) - corr_center_aZ2 )./ corr_radii_aZ2;

mo(:,2:4) = acc1(:,1:3);
mo(:,5:7) = acc2(:,1:3);
mo(:,11:13) = m(:,1:3);

ax=acc1(:,1);
ay=acc1(:,2);
az=acc1(:,3);
ax2=acc2(:,1);
ay2=acc2(:,2);
az2=acc2(:,3);
aGo=sqrt(mo(:,2).^2+mo(:,3).^2+mo(:,4).^2);
aG2o=sqrt(mo(:,5).^2+mo(:,6).^2+mo(:,7).^2);
aG=sqrt(ax.^2+ay.^2+az.^2);
aG2=sqrt(ax2.^2+ay2.^2+az2.^2);
%acc1、acc2、mag1 Calibrated

GG=sqrt(mo(:,8).^2+mo(:,9).^2+mo(:,10).^2)*0.01;

mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);
% mm42=sqrt(m2(:,1).^2+m2(:,2).^2+m2(:,3).^2);
% mm43=sqrt(m3(:,1).^2+m3(:,2).^2+m3(:,3).^2);
% mm44=sqrt(m4(:,1).^2+m4(:,2).^2+m4(:,3).^2);
% mm45=sqrt(m5(:,1).^2+m5(:,2).^2+m5(:,3).^2);
% mm46=sqrt(m6(:,1).^2+m6(:,2).^2+m6(:,3).^2);


x=smooth(m(:,1),smk);
y=smooth(m(:,2),smk);
z=smooth(m(:,3),smk);
smedm=[x y z];
xyz=sqrt(x.^2+y.^2+z.^2);
%% CHECK TIME
close all

check=[14600 18900 25000 35000];
[~,ck]=size(check);
ct=[];
for i=1:ck
ct(end+1)=t(check(i)); % transfer to real time
end

ticks=[check];
ticks=reshape (ticks, 1, numel(ticks));
ticks=sort(ticks,'ascend');

hold on
plot(m)
plot(aG)
set(gca, 'xtick', ticks);
grid on

%% power
close all
powerdata=xlsread('proj/7 2 中油 iron dic and pull/POWER3.CSV');
powerandtime=[(powerdata(:,1)/1000) powerdata(:,3)];
delemagn='定點8';

% 須找到mt0點時間
% mt0=509.5; %7 2 中油 power2.csv %%膠布脫落%%
mt0=1967; %7 2 中油 power3.csv 1%LOG46
% mt0=4273; %7 2 中油 power4.csv 2
% mt0=790; %7 2 中油 power6.csv %LOG51 3(5.4m power才有記錄)
% mt0=2453.5; %7 2 中油 power7.csv 4 deletp=30
% mt0=4500.5; %7 2 中油 power8.csv 5
% mt0=6167.2; %7 2 中油 power9.csv 6
% mt0=7724.5; %7 2 中油 power10.csv %LOG52 fix 7
% mt0=9294.5; %7 2 中油 power11.csv 8

onpt=0;
offpt=[];

[pk,~]=size(powerandtime);
for i=1:pk
    if powerandtime(i,2)==1 
        onpt(end+1)=powerandtime(i,1);
    end
    if powerandtime(i,2)==0
        offpt(end+1)=powerandtime(i,1);
    end
end
offpt(end+1)=9^9;

deletp=0;
onpt(1:deletp)=[];
offpt(1:deletp)=[];

onpt=onpt+mt0;
offpt=offpt+mt0;

ptstart=[];
ptend=[];
[~,onptk]=size(onpt);
for i=2:onptk
    if (onpt(i)-onpt(i-1))>2.1 % 區別出每段數據
        ptstart(end+1)=onpt(i);  
    end
    if (offpt(i)-offpt(i-1))>2.1  
        ptend(end+1)=offpt(i-1);
    end
end

tonpt=[];
toffpt=[];
[~,tonptk]=size(onpt);
for i=1:tonptk
[~,tonpt(i)]=min(abs(t-onpt(i)));
[~,toffpt(i)]=min(abs(t-offpt(i)));
end

tptstart=[];
tptend=[];
[~,ptstartk]=size(ptstart);
for i=1:ptstartk
[~,tptstart(i)]=min(abs(t-ptstart(i)));
[~,tptend(i)]=min(abs(t-ptend(i)));
end
extrat=0;
tptstart=tptstart-extrat*100;
tptend=tptend+extrat*100;
%%
tptstart=[43800 50300 54750 59850 72050 79000 87000];
tptstart=tptstart-600;
tptend=tptstart+1200;
[~,ptstartk]=size(tptstart);


% ticks=[tonpt toffpt];
ticks=[tptstart tptend];
ticks=reshape (ticks, 1, numel(ticks));
ticks=sort(ticks,'ascend');

hold on
% plot(mm4)
plot(m)
% plot(aG)
set(gca, 'xtick', ticks);
grid on

x=smooth(m(:,1),smk);
y=smooth(m(:,2),smk);
z=smooth(m(:,3),smk);
xyz=sqrt(x.^2+y.^2+z.^2);

close all

waveN=3;
magmin(waveN,ptstartk-1)=0;
magmax(waveN,ptstartk-1)=0;

useddata=y; % 用來判斷的數據

magmin=[];
magmax=[];
hold on
plot((aG(t0(1):t0(2))-0.5))
% plot(useddata)
plot(m)
filtCutOff = 0.8; % 高通判斷數據使可判斷波峰谷
[b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
useddataL = filtfilt(b, a, useddata);
plot(useddataL)
for i=1:ptstartk %找出每段的波峰波谷
testdata=useddataL(tptstart(i):tptend(i));
IndMin=find(diff(sign(diff(testdata)))>0)+20;  %獲得局部最小值的位置(可微調位置
IndMax=find(diff(sign(diff(testdata)))<0)+20;  %獲得局部最大值的位置

IndMin(waveN+1:size(IndMin))=[]; % 使同長度
IndMax(waveN+1:size(IndMax))=[];

magmin(:,i)=IndMin+tptstart(i); % 彌補時間差
magmax(:,i)=IndMax+tptstart(i);

plot(magmin(:,i),useddata(magmin(:,i)),'r^') % 標示在原數據
plot(magmax(:,i),useddata(magmax(:,i)),'k*')
end

% driftk=70; %7 2 中油 定點1
% magmax(:,17)=magmax(:,17)+driftk; % 飄移用
% magmin(:,17)=magmin(:,17)+driftk;

tag=[magmax magmin]; % 檢查位置
tag=reshape (tag, 1, numel(tag));
tag=sort(tag,'ascend');
set(gca, 'xtick', tag);
grid on
hold off
xlabel('Time')
ylabel('mag、acc')


ew(ptstartk,:)=0;
for i=1:ptstartk %每個位置磁力平均值
    for j=1:waveN
        wx(j,i)=x(magmax(j,i))-x(magmin(j,i));
        wy(j,i)=y(magmax(j,i))-y(magmin(j,i));
        wz(j,i)=z(magmax(j,i))-z(magmin(j,i));
    end   
    
    ew(i)=sqrt(mean(wx(:,i))^2+mean(wy(:,i))^2+mean(wz(:,i))^2);
    ews(i,:)=[mean(wx(:,i)) mean(wy(:,i)) mean(wz(:,i))];
    
%         if(i>=10&&i<=17) %定1
%         if(i>=8&&i<=18) %定2
%         if(i>=1&&i<=12) %定3
%         if(i>=5&&i<=22) %定4
%         if(i>=4&&i<=23) %定5
%         if(i>=2&&i<=25) %定6
%         if(i>=1&&i<=27) %定7
%         if(i>=1&&i<=29) %定8
% 
%          ews(i,:)=ews(i,:)*-1;
%          ew(i,:)=ew(i,:)*-1;
%         end
end
ew
%%
onoffacc1(ptstartk,3)=0; 
onoffacc2(ptstartk,3)=0; 
for i=1:ptstartk  %每個位置加速度
        onoffacc1(i,:)=mean(acc1(ptstart(i):ptend(i),:));
        onoffacc2(i,:)=mean(acc2(ptstart(i):ptend(i),:));
end
%%
delemagdata=[[1:i]' ews onoffacc1 onoffacc2 ];
text={'position' 'mx' 'my' 'mz' 'a1x' 'a1y' 'a1z' 'a2x' 'a2y' 'a2z'};
xlswrite('test1.xlsx',text,delemagn)
xlswrite('test1.xlsx',delemagdata,delemagn,'A2')
%% fixm(單點)
% fixm=[211008;211108]; %7 2 中油 定點1
fixm=[458746 458945;458647 458846]; %7 2 中油 定點2

fwx=[];fwy=[];fwz=[];
[~,fk]=size(fixm);
for i=1:fk
        fwx(i)=x(fixm(1,i))-x(fixm(2,i));
        fwy(i)=y(fixm(1,i))-y(fixm(2,i));
        fwz(i)=z(fixm(1,i))-z(fixm(2,i));
end
facc1=mean(acc1(fixm(1,1):fixm(2,fk),:));
facc2=mean(acc2(fixm(1,1):fixm(2,fk),:));
fix=[mean(fwx) mean(fwy) mean(fwz) facc1 facc2]
%% pull
pt0=0;
% pt0=668.3400; %5 5iron
pt0=205; %7 2 中油 blank 1、8%LOG53
pt0=1110; %7 2 8 pull
pt0=1311; %7 2 8 blank 2
pt0=1430; %7 2 8 pull 2
pt0=1699; %7 2 8 pull 3
pt0=2032; %7 2 9 blank 1
pt0=2195; %7 2 9 pull 1
pt0=137; %7 2 9 pull 2 %LOG55
% pt0=290; %7 2 9 pull 3
pont=[];pofft=[];
% ponmin=[48 56 61 64 69 72];% 5 5 iron pull
% ponsec=[51 52 42 52 13 20];
% poffmin=[50 58 63 67 70 74];
% poffsec=[38 25 21 23 2 5];
% pont=[0;-1];%7 2 中油 blank 1 %LOG53
% pofft=[0;47];
% pont=[0;-1];%7 2 中油 blank 8
% pofft=[0;50];
% pont=[0;0];%7 2 中油 8 pull
% pofft=[0;30];
% pont=[0;0];%7 2 中油 8 blank 2
% pofft=[0;20];
% pont=[0;0];%7 2 中油 8 pull 2
% pofft=[0;20];
% pont=[0;0];%7 2 中油 8 pull 3
% pofft=[0;20];
% pont=[0;0];%7 2 中油 8 pull 4
% pofft=[0;15];
% pont=[0;0];%7 2 中油 8 pull 5
% pofft=[0;20];
% pont=[0;-1];%7 2 中油 9 pull 2 %LOG55
% pofft=[0;20];
pont=[0;-1];%7 2 中油 9 pull 3
pofft=[0;20];


[~,k]=size(pont);
for i=1:k
pon(i)=pont(1,i)*60+pont(2,i)+pt0; %%%%%%%%%%%
poff(i)=pofft(1,i)*60+pofft(2,i)+pt0;
end
[~,k]=size(pon);
pton(k,:)=0;
ptoff(k,:)=0;
for i=1:k
[~,pton(i)]=min(abs(t-pon(i)));
[~,ptoff(i)]=min(abs(t-poff(i)));
end
ticks=[];
ticks=[pton ptoff];
ticks=reshape (ticks, 1, numel(ticks));
ticks=sort(ticks,'ascend');
hold on
plot(m)
% plot(aGsq)
% plot(aG*0.1)
set(gca, 'xtick', ticks);
grid on
xlabel('Time')
ylabel('Gauss') 
%% 拉動 data
filename='點 8 pull 7';
for i=1:1
pulln=int2str(i);
% pulldata=[t(pton(i):ptoff(i)) mo((pton(i):ptoff(i)),2:4) mo((pton(i):ptoff(i)),5:7) mo((pton(i):ptoff(i)),8:10) mo((pton(i):ptoff(i)),11:13) endcoder(pton(i):ptoff(i))-endcoder(pton(i)) h1(pton(i):ptoff(i)) h2(pton(i):ptoff(i))];
% text={'time' 'a1x' 'a1y' 'a1z' 'a2x' 'a2y' 'a2z' 'gx' 'gy' 'gz' 'mx' 'my' 'mz' 'encoder' 'encodero1' 'encodero2'};
pulldata=[t(pton(i):ptoff(i)) mo((pton(i):ptoff(i)),2:4) mo((pton(i):ptoff(i)),5:7) mo((pton(i):ptoff(i)),8:10) mo((pton(i):ptoff(i)),11:13)];
text={'time' 'a1x' 'a1y' 'a1z' 'a2x' 'a2y' 'a2z' 'gx' 'gy' 'gz' 'mx' 'my' 'mz'};
xlswrite(filename,text,pulln)
xlswrite(filename,pulldata,pulln,'A2')
end
%% 取得拉動平均速度
smf=300; %調整smooth參數判斷 
aGsq=abs(aG-mode(aG(:)));
aGsq=smooth(aGsq(t0(1):t0(2)),smf);
stationary=aGsq>0.01; %設定閥值 使0為靜止 1為移動 
close  all
hold on
% plot(stationary)
plot(aGsq)
plot(mm4)

for i=1:k   
    at(i,1)=find(stationary(pton(i):ptoff(i)) > 0, 1, 'first')+pton(i);
    at(i,2)=find(stationary((pton(i)+ptoff(i))/2:ptoff(i)) == 0, 1, 'first')+(pton(i)+ptoff(i))/2;
end
at=round(at)

atticks=[at(:,1) at(:,2) pton ptoff];
atticks=reshape (atticks, 1, numel(atticks));
atticks=sort(atticks,'ascend');
set(gca, 'xtick', atticks);
grid on
hold off

D=9.8; %已知移動距離(m)
pv=zeros(k,1);
 for i=1:k %  距離/移動時間
     pv(i,1)=D/(t(at(i,2))-t(at(i,1)))
 end
   
%% hall %%下多似計步器推算
addpath('hall');
H=mo(t0(1):t0(2),:); %因以加速度濾波及編碼器計算 須修正至量測區間
[Nf,~]=size(H);
h1=H(:,29);
h2=H(:,30);
% h1=smooth(h1,10);
% h2=smooth(h2,10);

map=[-10 10]; %map係數
mh1=(h1-min(h1))*(map(2)-map(1))/(max(h1)-min(h1))+map(1);
mh2=(h2-min(h2))*(map(2)-map(1))/(max(h2)-min(h2))+map(1);

close all
hold on
plot(mh1)
plot(mh2)
plot(mm4(t0(1):t0(2),:)*100)
plot(aG(t0(1):t0(2),:))

t=H(:,1);
tt=t-t(1);

inputEncoder1 = [tt mh1];
inputEncoder2 = [tt mh2];
Encoder_threshold=10; % 參數須設定

close all
[Encoder_X1,turning_velocity1,turning_degree1] = Encoder_1(inputEncoder1,Encoder_threshold);

close all
[Encoder_X2,turning_velocity2,turning_degree2] = Encoder_1(inputEncoder2,Encoder_threshold);
Encoder_X1(Nf)
Encoder_X2(Nf)
endcoder=Encoder_X1;
if endcoder(Nf)<Encoder_X2(Nf)
    endcoder=Encoder_X2;
end    
%%
close  all

figure(1)
hold on
plot(inputEncoder1(:,1),Encoder_X1);
plot(inputEncoder2(:,1),Encoder_X2);
legend('hall1','hall2','Location','southeast')
title('Encoder')
hold off

% figure(2)
% hold on
% plot(endcoder,smooth(H(:,11),25)); %m的x軸
% % plot(endcoder,smooth(m6(:,1),25));
% xlabel('postion(m) (hall2)')
% ylabel('magnetic field strength(guess) (x)')
% legend('mag1','mag6','Location','southeast')
% title('magnetic field strength in each position')
% hold off

% stopt=[1 200];  %% 設定靜止時間(資料數)
% dt=0.01;
% [M_velo2, dis_from_Mag2]=getVtromMag(m6*100, m*100, dt, N ,stopt);
% M_dis2=cumtrapz(tt,M_velo2);
% M_dis2(N)
% 
% 
% figure(3) %hall速度 vs 磁場推算速度
% hold on
% plot(tt,Encoder_X2);
% plot(tt,M_dis2);
% xlabel('time(s)')
% ylabel('postion(m)  (hall2)')
% legend('hall2 sensor','getVfromMag','Location','southeast')
% title('Encoder position and getVfromMag position')
% hold off
%                acc1      acc2      gy       mag
fuckyoutwo=[tt H(:,2:4) H(:,5:7) H(:,8:10) H(:,11:13) endcoder h1 h2];
text={'time' 'a1x' 'a1y' 'a1z' 'a2x' 'a2y' 'a2z' 'gx' 'gy' 'gz' 'mx' 'my' 'mz' 'encoder' 'encodero1' 'encodero2'};
% xlswrite('test.xlsx',text,testn)
% xlswrite('test.xlsx',fuckyoutwo,testn,'A2')
% sprintf(testn)
%% 編碼器確認位移
for i=2:k
    encd(i)=Encoder_X1(tect(i))-Encoder_X1(tect(i-1));
    encd2(i)=Encoder_X2(tect(i))-Encoder_X2(tect(i-1));
end
encd'
encd2'
hold on
plot(Encoder_X1)
plot(Encoder_X2)
set(gca, 'xtick', ticks(24));
grid on
%% threshold filter ((fixing
smf=265; %調整smooth參數判斷 
aGsq=abs(aG-mode(aG(:)));
aGsq=smooth(aGsq(t0(1):t0(2)),smf);
stationary=aGsq>0.0178; %0為靜止 1為移動 
close  all
hold on
plot(stationary)
plot(aGsq)
hold off
%%

T0=[];T1=[];
c=0; %condition
fixt=50; %微調區間
x=smooth(H(:,11),smk);
y=smooth(H(:,12),smk);
z=smooth(H(:,13),smk);
xyz=sqrt(x.^2+y.^2+z.^2);

for i=2:Nf %判斷移動靜止區間
    if stationary(i)==1 && c==0
        T0(end+1)=i;
        c=1;
    end
    if stationary(i)==0 && c==1
        T1(end+1)=i+fixt;
        c=2;
    end
    if stationary(i)==1 && c==2
        T0(end+1)=i;
        c=1;
    end
end
T1=[T0(1)-500 T1(:,1:end)];
T0=[T0(:,1:end) T1(end)+500];

tticks=[T0 T1];
% tticks=[T0];
tticks=reshape (tticks, 1, numel(tticks));
tticks=sort(tticks,'ascend');

close all
hold on
plot(aG(t0(1):t0(2))-0.5)
plot(aGsq)
plot(xyz)
% plot(xyz2)
% plot(y)
set(gca, 'xtick', tticks);
grid on
%%
acc=H(:,2:4);
samplePeriod=1/100;
acc = acc * 9.81;

% Integrate acceleration to yield velocity 整合加速度以產生速度
vel = zeros(size(acc));
for t = 2:length(vel)
    vel(t,:) = vel(t-1,:) + acc(t,:) * samplePeriod;
    if(stationary(t) == 0)
        vel(t,:) = [0 0 0];     % force zero velocity when foot stationary 當腳靜止時強制零速度
    end
end
pos = zeros(size(vel));
for t = 2:length(pos)
    pos(t,:) = pos(t-1,:) + vel(t,:) * samplePeriod;    % integrate velocity to yield position 將速度整合到屈服位置
end
poss=sqrt(pos(:,1).^2+pos(:,2).^2+pos(:,3).^2)*0.01;
poss(t)
plot(poss)


%% threshold filter waves (for pipe measuring)
close all
[~,k]=size(T0);

waveN=1;
magmin(waveN,k-1)=0;
magmax(waveN,k-1)=0;

useddata=y; % 用來判斷的數據

hold on
plot((aG(t0(1):t0(2))-0.5))
plot(useddata)
filtCutOff = 1; % 高通判斷數據使可判斷波峰谷
[b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
useddataL = filtfilt(b, a, useddata);
plot(useddataL)
for i=1:k %找出每段的波峰波谷
testdata=useddataL(T1(i):T0(i));
[~,kt]=size(testdata);
IndMin=find(diff(sign(diff(testdata)))>0)+20;  %獲得局部最小值的位置(可微調位置
IndMax=find(diff(sign(diff(testdata)))<0)+20;  %獲得局部最大值的位置

IndMin(waveN+1:size(IndMin))=[]; % 使同長度
IndMax(waveN+1:size(IndMax))=[];

magmin(:,i)=IndMin+T1(i); % 彌補時間差
magmax(:,i)=IndMax+T1(i);

plot(magmin(:,i),useddata(magmin(:,i)),'r^') % 標示在原數據
plot(magmax(:,i),useddata(magmax(:,i)),'k*')
end
plot(endcoder)

tag=[T0  T1]; % 取樣區間
% tag=[magmax(1,:) magmin(1,:)]; % 檢查位置
tag=reshape (tag, 1, numel(tag));
tag=sort(tag,'ascend');
set(gca, 'xtick', tag);
grid on
hold off
xlabel('Time')
ylabel('mag、acc、encoder')


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

waveEncoder(k,:)=0; 
for i=1:k  %每個位置編碼器的平均值
        waveEncoder(i)=mean(endcoder(T1(i):T0(i)));
end
waveEncoder;
%% %%%%%% elecmag d %%%%%%
% figure(2)
% tt=mo(:,38); %frequency
%  plot(tt)
t=mo(:,1);
x=smooth(m(:,1),smk);y=smooth(m(:,2),smk);z=smooth(m(:,3),smk);
xyz=sqrt(x.^2+y.^2+z.^2);
%%%%%
mins=[0 5 7 7 8 9 11 12 13 14 15 16 17 18 19 20 22 24 24 25 27 28 29 30 31 32 33 34 35 36 38 39]';% 5 11 對稱
secs=[30 30 0 50 40 55 10 25 10 10 10 20 20 20 20 0 20 0 40 50 0 0 0 1 0 0 10 30 50 50 20 40]';

[k,~]=size(mins);
for i=1:k
onn(i)=mins(i)*60+secs(i);
end

on=[t(13500) t(16700) t(18800) t(21350) t(23400) t(26500) t(36000) t(42300)]; %air % 8 4 end huge wire
% on=[t(14400) t(16900) t(21000) t(23400) t(25400) t(29000) t(32800) t(35800)]; %iron
% on=[t(7300) t(13400) t(16700)]; %42v
% on=[t(14600) t(18900) t(25000) t(35000)]; % d high

% off=on+3
off=on-4;;

% on=1;
% off=2;

% ticks=[on off];
% ticks=sort(ticks,'ascend');
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

for i=1:k % use time
    ed(i)=sqrt((x(ton(i))-x(toff(i)))^2+(y(ton(i))-y(toff(i)))^2+(z(ton(i))-z(toff(i)))^2);
    eds(i,:)=[x(ton(i))-x(toff(i)) y(ton(i))-y(toff(i)) z(ton(i))-z(toff(i))]';
end
% for i=1:k % use dataN
%     ed(i)=sqrt((x(on(i))-x(off(i)))^2+(y(on(i))-y(off(i)))^2+(z(on(i))-z(off(i)))^2);
%     eds(i,:)=[x(on(i))-x(off(i)) y(on(i))-y(off(i)) z(on(i))-z(off(i))]';
% end

ed
% eds

% onoffEncoder(k,:)=0; 
% for i=1:k  %每個位置編碼器
%         onoffEncoder(i)=mean(endcoder(toff(i)-t0(1):ton(i)-t0(1)));
% end
onoffacc1(k,3)=0; 
onoffacc2(k,3)=0; 
for i=1:k  %每個位置加速度
        onoffacc1(i,:)=mean(acc1(toff(i):ton(i),:));
        onoffacc2(i,:)=mean(acc2(toff(i):ton(i),:));
end

close all
% figure(1);
hold on
% % plot(d,ed)
% % plot(y)
% % plot(m(:,1))
% % plot(m(:,2))
% % plot(m(:,3))
plot(x)
plot(y)
plot(z)
% plot(mm4)
plot(xyz)

% plot(G)
% plot(A1)
plot(aG2)

tickson=[ton];
ticksoff=[toff];
tickson=reshape (tickson, 1, numel(tickson));
ticksoff=reshape (ticksoff, 1, numel(ticksoff));
tickson=sort(tickson,'ascend');
ticksoff=sort(ticksoff,'ascend');

for i=1:k
    plot(tickson(i),xyz(tickson(i)),'bo')
    plot(ticksoff(i),xyz(ticksoff(i)),'ro')
end
set(gca, 'xtick', ticks);
grid on
xlabel('Time')
ylabel('Gauss')
hold off
%% on off 定點data
delemagn='定點';
delemagdata=[[1:k]' onoffacc1 onoffacc2 eds];%onoffEncoder];
text={'position' 'a1x' 'a1y' 'a1z' 'a2x' 'a2y' 'a2z' 'mx' 'my' 'mz'};% 'encoder'};
xlswrite('test.xlsx',text,delemagn)
xlswrite('test.xlsx',delemagdata,delemagn,'A2')
%% iron wave 6magggg (toff)

tags(k,2)=0;
atag(k,:)=0;
atag1(k,:)=0;
waveN=4;
magmin(waveN,k)=0;
magmax(waveN,k)=0;

useddata=xyz; % 用來判斷的數據

hold on
plot(useddata)
filtCutOff = 0.4; % 高通判斷數據使可判斷波峰谷 %0.9
[b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
useddataL = filtfilt(b, a, useddata);
plot(useddataL)
for i=1:k
testdata=useddataL([ton(i):toff(i)],:);
% testdata=xyz([toff(i):ton(i)],:);

tag=ones(size(on))*100;
tagg=ones(size(tag))*1200;



tags(:,1)=tag;
tags(:,2)=tagg;

testdata(1:tags(i,1))=[];%
testdata(tags(i,2):size(testdata))=[];%
% figure(i);
% hold on
% plot(testdata)
[b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
testdata = filtfilt(b, a, testdata);
% plot(testdata)
c=findpeaks(testdata);
IndMin=find(diff(sign(diff(testdata)))>0)+22;%1.5off +28   %獲得局部最小值的位置
IndMax=find(diff(sign(diff(testdata)))<0)+22;   %獲得局部最大值的位置
% figure; hold on; box on;
% plot(1:length(testdata),testdata);
% plot(IndMin,testdata(IndMin),'r^')
% plot(IndMax,testdata(IndMax),'k*')
% legend('原數據','低通','波谷','波峰')

atag(i)=ton(i)+tags(i,1);
atag1(i)=ton(i)+tags(i,1)+tags(i,2);

% if(i==2)
% magmin(:,i)=IndMin+atag(i);
% magmax(:,i)=IndMax+atag(i);
% plot(magmin(:,i),xyz(magmin(:,i)),'r^')
% plot(magmax(:,i),xyz(magmax(:,i)),'k*')
    
% else
    IndMin(waveN+1:size(IndMin))=[]; % 使同長度
    IndMax(waveN+1:size(IndMax))=[];
% end

magmin(:,i)=IndMin+atag(i);
magmax(:,i)=IndMax+atag(i);
plot(magmin(:,i),useddata(magmin(:,i)),'r^')
plot(magmax(:,i),useddata(magmax(:,i)),'k*')
end

plot(toff,useddata(toff),'bo')

mt=m; %用來計算結果的數據
xt=smooth(mt(:,1),smk);yt=smooth(mt(:,2),smk);zt=smooth(mt(:,3),smk);
xyzt=sqrt(xt.^2+yt.^2+zt.^2);
plot(xt)
plot(yt)
plot(zt)
% plot(xyzt)
plot(G)
% plot(A1)
% plot(A2)

% tag=[atag atag1];
wtag=[magmin magmax];
% tag=reshape (tag, 1, numel(tag));
wtag=reshape (wtag, 1, numel(wtag));
% tag=[tag wtag ticks];
tag=[wtag ticks];
tag=reshape (tag, 1, numel(tag));
tag=sort(tag,'ascend');
set(gca, 'xtick', tag);
grid on
xlabel('Time')
ylabel('Gauss')
legend('原數據','低通','波谷','波峰')

hold off

ew(k,:)=0;
ews(k,3)=0;
wxb(waveN,k)=0;
wyb(waveN,k)=0;
wzb(waveN,k)=0;
for i=1:k % wave ,use time
%     for j=1:waveN %波峰減toff
% %         wxb(j,i)=xt(magmax(j,i));
% %         wyb(j,i)=yt(magmax(j,i));
% %         wzb(j,i)=zt(magmax(j,i));
%         wxb(j,i)=xt(magmin(j,i));
%         wyb(j,i)=yt(magmin(j,i));
%         wzb(j,i)=zt(magmin(j,i));
%     end
%     ew(i)=sqrt((mean(wxb(:,i))-xt(toff(i)))^2+(mean(wyb(:,i))-yt(toff(i)))^2+(mean(wzb(:,i))-zt(toff(i)))^2);
%     ews(i,:)=[mean(wxb(:,i))-xt(toff(i)) mean(wyb(:,i))-yt(toff(i)) mean(wzb(:,i))-zt(toff(i))];
% % %     if(i==3|i==5|i==6|i==7|i==8|i==9|i==15|i==16)%max 1224 201
% % %     if(i~=3&i~=5&i~=6&i~=7&i~=8&i~=9&i~=15&i~=16)%min
% % %     if(i==6|i==7|i==8|i==9|i==16)%max 1224 160.5
% % %     if(i~=6&i~=7&i~=8&i~=9&i~=16)%min
% % %     if(i==8|i==9|i==10|i==17)%max 1224 120
% % %     if(i~=8&i~=9&i~=10&i~=17)%min
% % %     if(i==7|i==8)%max 1225 pvc201
% % %     if(i~=7&i~=8)%min
% % % 
% % %          ews(i,:)=ews(i,:)*0;
% % %          ew(i,:)=ew(i,:)*0;
% % %     end



    for j=1:waveN  %波峰減波谷
        wx(j,i)=xt(magmax(j,i))-xt(magmin(j,i));
        wy(j,i)=yt(magmax(j,i))-yt(magmin(j,i));
        wz(j,i)=zt(magmax(j,i))-zt(magmin(j,i));
    end
    ew(i)=sqrt(mean(wx(:,i))^2+mean(wy(:,i))^2+mean(wz(:,i))^2);
    ews(i,:)=[mean(wx(:,i)) mean(wy(:,i)) mean(wz(:,i))];
end
ew
ews
%%
for i=1:k
    toffm(i,:)=[xt(toff(i)) yt(toff(i)) zt(toff(i))];
end
toffm

%% wave magggg bro

tags(k,2)=0;
atag(k,:)=0;
atag1(k,:)=0;
waveN=8;
magmin(waveN,k)=0;
magmax(waveN,k)=0;


hold on
plot(xyz)
filtCutOff = 0.4;%0.9
[b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
xyzL = filtfilt(b, a, xyz);
plot(xyzL)
for i=1:k
testdata=xyz([ton(i):toff(i)],:);
% testdata=xyz([toff(i):ton(i)],:);
% tag=[250 200 270 160 140];%11/14
% tag=[50 50 50 50 50 50];%
% tag=[100 200 200 200 100 100 200 200 200 200];%
% tag=[100 200 200 200 100 200 200 200];%

% tag=[100 300 150 300 150 450 350 350 350 350 450 450 550 150 150];%11/19

% tag=ones(size(on))*600;
% tagg=ones(size(tag))*3500;
% tag=ones(size(on))*50;
% tagg=ones(size(tag))*3300;
tag=ones(size(on))*300; % 11 28 iron
tagg=ones(size(tag))*3300;
% tag=ones(size(on))*150; % 11 29 2elec mag
% tagg=ones(size(tag))*4400;
% tag=ones(size(on))*200; %iron 12 9 1M(no)
% tagg=ones(size(tag))*3300;

% tagg=ones(size(tag))*3000; % 11 26 (on*1)

tags(:,1)=tag;
tags(:,2)=tagg;

testdata(1:tags(i,1))=[];%
testdata(tags(i,2):size(testdata))=[];%
% figure(i);
% hold on
% plot(testdata)
[b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
testdata = filtfilt(b, a, testdata);
% plot(testdata)
c=findpeaks(testdata);
IndMin=find(diff(sign(diff(testdata)))>0)+22;   %獲得局部最小值的位置
IndMax=find(diff(sign(diff(testdata)))<0)+22;   %獲得局部最大值的位置
% figure; hold on; box on;
% plot(1:length(testdata),testdata);
% plot(IndMin,testdata(IndMin),'r^')
% plot(IndMax,testdata(IndMax),'k*')
% legend('原數據','低通','波谷','波峰')

atag(i)=ton(i)+tags(i,1);
atag1(i)=ton(i)+tags(i,1)+tags(i,2);
% atag(i)=toff(i)+tags(i,1);
% atag1(i)=toff(i)+tags(i,1)+tags(i,2);

IndMin(waveN+1:size(IndMin))=[]; % 使同長度
IndMax(waveN+1:size(IndMax))=[];

magmin(:,i)=IndMin+atag(i);
magmax(:,i)=IndMax+atag(i);
plot(magmin(:,i),xyz(magmin(:,i)),'r^')
plot(magmax(:,i),xyz(magmax(:,i)),'k*')
end
    
plot(x)
plot(y)
plot(z)
% plot(G)
% plot(A1)
% plot(A2)

tag=[atag atag1];
tag=reshape (tag, 1, numel(tag));
tag=sort(tag,'ascend');
set(gca, 'xtick', tag);
grid on
hold off

set(gca, 'xtick', ticks);
grid on
xlabel('Time')
ylabel('Gauss')
% legend('原數據','低通','波谷','波峰','','','','','','','','','','','','','','','x','y','z')

ew(k,:)=0;
ews(k,3)=0;
wx(waveN,k)=0;
wy(waveN,k)=0;
wz(waveN,k)=0;
for i=1:k % wave ,use time
%     for j=1:waveN
%         wx(j,i)=x(magmin(j,i))-x(magmax(j,i));
%         wy(j,i)=y(magmin(j,i))-y(magmax(j,i));
%         wz(j,i)=z(magmin(j,i))-z(magmax(j,i));
%     end
%     ew(i)=sqrt(mean(wx(:,i))^2+mean(wy(:,i))^2+mean(wz(:,i))^2);
%     ews(i,:)=[mean(wx(:,i)) mean(wy(:,i)) mean(wz(:,i))];
% %     if(i==4|i==7|i==8|i==9|i==10|i==11|i==12|i==13)
%     if(i>0)
% 
%         ews(i,:)=ews(i,:)*-1;
%     end
    
        for j=1:waveN   %toff
        wxb(j,i)=xt(magmax(j,i));
        wyb(j,i)=yt(magmax(j,i));
        wzb(j,i)=zt(magmax(j,i));
%         wxb(j,i)=xt(magmin(j,i));
%         wyb(j,i)=yt(magmin(j,i));
%         wzb(j,i)=zt(magmin(j,i));
    end
    ew(i)=sqrt((mean(wxb(:,i))-xt(toff(i)))^2+(mean(wyb(:,i))-yt(toff(i)))^2+(mean(wzb(:,i))-zt(toff(i)))^2);
    ews(i,:)=[mean(wxb(:,i))-xt(toff(i)) mean(wyb(:,i))-yt(toff(i)) mean(wzb(:,i))-zt(toff(i))];

end
ew
ews
%% 1 wave mag test
% hi=1;
% testdata=xyz([ton(hi):toff(hi)],:);
% tag=[1 600];
% testdata(1:tag(1))=[];
% testdata(tag(2):size(testdata))=[];
% % figure(i);
% hold on
% plot(testdata)
% filtCutOff = 2;
% [b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
% testdata = filtfilt(b, a, testdata);
% plot(testdata)
% c=findpeaks(testdata);
% IndMin=find(diff(sign(diff(testdata)))>0)+1;   %?得局部最小值的位置
% IndMax=find(diff(sign(diff(testdata)))<0)+1;   %?得局部最大值的位置
% % figure; hold on; box on;
% plot(1:length(testdata),testdata);
% plot(IndMin,testdata(IndMin),'r^')
% plot(IndMax,testdata(IndMax),'k*')
% % legend('原數據','低通','波谷?','波峰?')

% tags=[];
% tags(2,:)=0;
% tags(1)=ton(hi)+tag(1);
% tags(2)=ton(hi)+tag(1)+tag(2);
% hold on
% plot(xyz)
% [b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
% xyzL = filtfilt(b, a, xyz);
% plot(xyzL)
% 
% IndMin=IndMin+tag(1)+ton(hi);
% IndMax=IndMax+tag(1)+ton(hi);
% plot(IndMin,xyz(IndMin),'r^')
% plot(IndMax,xyz(IndMax),'k*')
% tags=reshape (tags, 1, numel(tags));
% tags=sort(tags,'ascend');
% set(gca, 'xtick', tags);
% grid on
%% 2 mag waves
% waveAb(waveN/2,k)=0;
% waveAs(waveN/2,k)=0;
% waveBb(waveN/2,k)=0;
% waveBs(waveN/2,k)=0;
% WS=[];
% ewbsA(k,:)=0;
% ewbssA(k,3)=0;
% ewbsB(k,:)=0;
% ewbssB(k,3)=0;
% wxA(waveN/2,k)=0;
% wyA(waveN/2,k)=0;
% wzA(waveN/2,k)=0;
% wxB(waveN/2,k)=0;
% wyB(waveN/2,k)=0;
% wzB(waveN/2,k)=0;
% 
% for i=1:k
%     waveAb(i,:)=magmax(1,:)+(i-1)*500;
%     waveAs(i,:)= waveAb(i,:)+145;
%     waveBb(i,:)=waveAb(i,:)+250;
%     waveBs(i,:)=waveAb(i,:)+385;
% end
%     waveAb(:,8)=waveAb(:,8)-130; % FIX
%     waveAs(:,8)=waveAs(:,8)-130;    
%     waveBb(:,8)=waveBb(:,8)-130; % FIX
%     waveBs(:,8)=waveBs(:,8)-130;  
%     
% figure(9)
% hold on
% plot(xyz)
% plot(xyzL)
% for i=1:k
% plot(waveAb(:,i),xyz(waveAb(:,i)),'r*')
% plot(waveAs(:,i),xyz(waveAs(:,i)),'b^')
% plot(waveBb(:,i),xyz(waveBb(:,i)),'g*')
% plot(waveBs(:,i),xyz(waveBs(:,i)),'y^')
% end
% 
% WS=[waveAb waveAs waveBb waveBs];
% WS=reshape (WS, 1, numel(WS));
% WS=sort(WS,'ascend');
% 
% legend('數據','數據','ab','as','bb','bs')
% set(gca, 'xtick', WS);
% grid on
% hold off
% 
% for i=1:k % wave ,use time
%     for j=1:waveN/2
%         wxA(j,i)=x(waveAb(j,i))-x(waveAs(j,i));
%         wyA(j,i)=y(waveAb(j,i))-y(waveAs(j,i));
%         wzA(j,i)=z(waveAb(j,i))-z(waveAs(j,i));
%     end
%     for j=1:waveN/2
%         wxB(j,i)=x(waveBb(j,i))-x(waveBs(j,i));
%         wyB(j,i)=y(waveBb(j,i))-y(waveBs(j,i));
%         wzB(j,i)=z(waveBb(j,i))-z(waveBs(j,i));
%     end
%     ewbsA(i)=sqrt(mean(wxA(:,i))^2+mean(wyA(:,i))^2+mean(wzA(:,i))^2);
%     ewbssA(i,:)=[mean(wxA(:,i)) mean(wyA(:,i)) mean(wzA(:,i))];
%     ewbsB(i)=sqrt(mean(wxB(:,i))^2+mean(wyB(:,i))^2+mean(wzB(:,i))^2);
%     ewbssB(i,:)=[mean(wxB(:,i)) mean(wyB(:,i)) mean(wzB(:,i))];
% %     if(i==4|i==7|i==8|i==9|i==10|i==11|i==12|i==13)
% %     if(i>0)
% 
% %         ewbss(i,:)=ewbss(i,:)*-1;
% %     end
% end
% ewbsA
% ewbsB
% % ewbssA
% % ewbssB
% 
% ewbsall=[ewbsB ewbsA];  %% wrong !
% ewbsall=reshape (ewbsall, 1, numel(ewbsall));
% % ewbsall=sort(ewbsall,'descend');
% % dws=[252.9822128 178.8854382 113.137085 80 113.137085 178.8854382 252.9822128 329.84845 200 144.222051 120 144.222051 200 268.3281573 341.7601498 417.6122604];
% % dws=sort(dws,'descend');
% % 
% % figure(10)
% % plot(dws,ewbsall)
% 
% figure(10)
% dws=[329.84845 252.9822128 178.8854382 113.137085 80 113.137085 178.8854382 252.9822128 417.6122604 341.7601498 268.3281573 200 144.222051 120 144.222051 200];
% hold on
% for i=1:k
%     plot(ewbsA(i),dws(i+8),'o')
%     plot(ewbsB(i),dws(i),'x')
% end
% 
% hold off
% 
% figure(11)
% loglog(dws,ewbsall)

%%
% pvc1=mo([20000:48000],:);
% % pvc2=mo([34000:60000],:);
% % pvc3=mo([96000:125000],:);
% % pvc4=mo([170000:194000],:);
% % pvc5=mo([217000:239000],:);
% % % 
% % iron=mo([16000:49300],:);
% % iron2=mo([55000:75000],:);
% 
% % pvc3empull=mo([36000:89000],:);
% % pvc3empull2=mo([93000:138000],:);
% 
% plot(mm4)
% 
% ttt=[20000 48000];
% % ttt=[34000 60000 96000 125000 170000 194000 217000 239000];
% % ttt=[16000 49300 55000 75000];
% % ttt=[36000 89000 93000 138000];
% 
% a1(k,3)=0;
% for i=1:k
%     a1(i,:)=sum(acc1([ton(i):toff(i)],:))/(abs(toff(i)-ton(i)))
% end
% 
% set(gca, 'xtick', ttt);
% grid on
%% threshold filter (for pipe measuring)(need fix

movet=500;
T0=[];T1=[];
c=0;
x=smooth(H(:,11),smk);y=smooth(H(:,12),smk);z=smooth(H(:,13),smk);
xyz=sqrt(x.^2+y.^2+z.^2);
% x2=smooth(m6(:,1),smk);y2=smooth(m6(:,2),smk);z2=smooth(m6(:,3),smk);
% xyz2=sqrt(x2.^2+y2.^2+z2.^2);

T0(end+1)=1;
T1(end+1)=t0(1)+300;
for i=t0(1):t0(2) %判斷移動後靜止的時間
    if aG(i)>aG(t0(1))*1.01 && c==0
        T0(end+1)=i-70;
        T1(end+1)=i+300;
        c=1;
    end
    if aG(i)>aG(t0(1))*1.01 && i>T0(end)+700
        T0(end+1)=i-70;
        T1(end+1)=i+300;
    end
end
T0(end+1)=t0(2)-200;
T1(end+1)=t0(2)-1;

tticks=[T0 T1];
% tticks=[T0];
tticks=reshape (tticks, 1, numel(tticks));
tticks=sort(tticks,'ascend');

close all
hold on
% plot(xyz)
% plot(xyz2)
plot(y)
plot(aG)
set(gca, 'xtick', tticks);
grid on

%% %%%%核康mag ? %%%%
% m2(1:361,:)=[]; %361 12v %2031 3.3v
% m2(1:8540,:)=[];m(1:8540,:)=[];
% m2(1:5421,:)=[];m(1:5421,:)=[];


% t0=400;t1=1000; %背景
% t2=2000;t3=2700; %開磁鐵區間
% x0=sum(m2([t0:t1],1))/(t1-t0); %取平均  D test
% y0=sum(m2([t0:t1],2))/(t1-t0);
% z0=sum(m2([t0:t1],3))/(t1-t0);
% x1=sum(m2([t2:t3],1))/(t3-t2);
% y1=sum(m2([t2:t3],2))/(t3-t2);
% z1=sum(m2([t2:t3],3))/(t3-t2);
% 
% x20=sum(m([t0:t1],1))/(t1-t0);%取平均  D test
% y20=sum(m([t0:t1],2))/(t1-t0);
% z20=sum(m([t0:t1],3))/(t1-t0);
% x21=sum(m([t2:t3],1))/(t3-t2);
% y21=sum(m([t2:t3],2))/(t3-t2);
% z21=sum(m([t2:t3],3))/(t3-t2);%
% 
% x2=x1-x0;
% y2=y1-y0;
% z2=z1-z0;
% x22=x21-x20;
% y22=y21-y20;
% z22=z21-z20;
% xyz0=[x0 y0 z0]';
% xyz1=[x1 y1 z1]';
% xyz2=[x2 y2 z2]';
% xyz20=[x20 y20 z20]';
% xyz21=[x21 y21 z21]';
% xyz22=[x22 y22 z22];
% m32=sqrt(x2^2+y2^2+z2^2);
% m33=sqrt(x22^2+y22^2+z22^2);
% liu=[x0 y0 z0 x1 y1 z1 x2 y2 z2 m32]';

% figure(1)

% hold on
% mm2=sqrt((m2(:,1)-x0).^2+(m2(:,2)-y0).^2+(m2(:,3)-z0).^2);
% mm22=sqrt((m(:,1)-x20).^2+(m(:,2)-y20).^2+(m(:,3)-z20).^2);
% mm42=sqrt(m2(:,1).^2+m2(:,2).^2+m2(:,3).^2);
% mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);
% plot(m2)
% plot(m)
% plot(mm22)
% plot(mm2)
% plot(mm4)
% plot(mm42)
% plot(m2(:,1))
% plot(m2(:,2))
% plot(m2(:,3))
% plot(m(:,1))
% plot(m(:,2))
% plot(m(:,3))
% scatter3(m2(:,1),m2(:,2),m2(:,3))
% scatter3(m(:,1),m(:,2),m(:,3))
% legend('一號','二號')
% hold off
%%

% sm0=std(mm2([t0:t1],1));
% sm1=std(mm2([t2:t3],1));
% sm=[sm0 sm1]';
% rotor=(mo2(:,17)+mo2(:,18))/2;
% t=mo2(:,1);
% lr=0;
% d(size(rotor))=0;
% nd(size(rotor))=0;
% v(size(rotor))=0;
% Rotor(size(rotor))=0;
% for i=1:size(rotor)
% r=rotor(i);
% Rotor(i)=Rotor(i)+rotor(i);
% if(lr==348.75&&r==0)
%     Rotor(1,(i:size(rotor)))=Rotor(1,(i:size(rotor)))+360;
% end
% lr=r;
% end
% Rotor=Rotor';
% r=4.52607796436989;
% 360*3980/2/pi/rotor(i);
% Rotor(i)/360*2*pi*r*0.01;
% rotor(i)/360*2*pi*r*0.01;
% for i=1:size(rotor)
%     d(i)=Rotor(i)/360*2*pi*r*0.01;
% end
% for i=2:size(rotor)
%     nd(i)=d(i)-d(i-1);
%     v(i)=(nd(i)-nd(i-1))/(t(i)-t(i-1));
% end
% ax=mo2(:,2)*9.8;
% ax=ax-(sum(ax([200:1200],1))/1000);
% 

%%
% a=[0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000]';
% b(size(a),:)=0;
% hold on
% data=xlload('elecmag d test\field test\data.xls');
% for j=1:9
% %     quiver(a,b,)
% end

% mson=[370973 333580 428659 504000];
% msoff=[364300 320000 421418 620000];
% % m300on=[333580];
% % m300off=[32000];
% % m160on=[370973];
% % m160off=[364300];
% % m500on=[428659];
% % m500off=[421418];
% ticks=[mson msoff];
% ticks=sort(ticks,'ascend');
% [~,k]=size(mson);
% ed(k,:)=0;
% for i=1:k
%     ed(i)=sqrt((x(mson(i))-x(msoff(i)))^2+(y(mson(i))-y(msoff(i)))^2+(z(mson(i))-z(msoff(i)))^2);
% end
% ed
% plot(mm4)
% plot(xyz)
% set(gca, 'xtick', ticks);
% grid on	
%% hall
% % H1=mo(:,14)/1000;
% % H2=mo(:,15)/1000;
% % H1=mo(:,14);
% % H2=mo(:,15);
% H1=smooth(mo(:,14)-500,1);
% H2=smooth(mo(:,15)-500,1);
% hold on
% 
% plot(H1)
% % filtCutOff = 0.5;
% % [b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
% % H1 = filtfilt(b, a, H1);
% % H2 = filtfilt(b, a, H2);
% 
% i=size(H1)-1;
% for i=2:i
% if(H1(i-1)==H1(i+1))&&(abs(H1(i)-H1(i-1))==1)&&(abs(H1(i)-H1(i+1))==1)
%     H1(i)=H1(i-1);
% i
% end
% 
% end
%     
% % plot(mm4)
% % xyz=xyz*10;
% % plot(xyz)
% % plot(x)
% % plot(y)
% % plot(z)
% plot(H1)
% % plot(H2)
%% iron test
% % oeon=[1.1 4.9 8.4 12 15.7 20 23.5 26.9 30.2 33.4]*1000; %%10 7
% oeon=[11.87 13.5 14.5 20.15 22.2]*1000; %%10 7 mother board
% % teon=[2.7 6.5 9.7 13.3 17 21.4 24.6 28.4 31.5 34.8]*1000;
% % omt=oeon-5*85;
% omt=oeon(1)-10*1000;
% 
% hold on
% plot(x)
% plot(y)
% plot(z)
% % plot(mm4)
% % plot(xyz)
% 
% tickss=[oeon omt];
% tickss=sort(tickss,'ascend');
% set(gca, 'xtick', tickss);
% grid on	
% xlabel('Time')
% ylabel('Gauss')
% 
% [~,k]=size(oeon);
% om(k,3)=0;
% % em(k,3)=0;
% em(k,:)=0;
% tm(k,3)=0;
% for i=1:k
% %     om(i,:)=[x(omt(i)) y(omt(i)) z(omt(i))];
% %     em(i,:)=[x(oeon(i)) y(oeon(i)) z(oeon(i))];
% %     tm(i,:)=[x(teon(i))-x(omt(i)) y(teon(i))-y(omt(i)) z(teon(i))-z(omt(i))];
% %     em(i)=sqrt((x(oeon(i))-x(omt))^2+(y(oeon(i))-y(omt))^2+(z(oeon(i))-z(omt))^2);
% end
% om
% em
% % tm


%%
%%%%%%%%%%% wire m2 d test%%%%%%%%%%
% d=[20 30 40 50 60 70 80 90 100 110 120 130 140 150]; %%  wire~~~~
% d=[230 220 210 200 190 180 170 160 150 140 130 120 110 100 90 80 70 60 50 40 30 20]; %%  wire~~~~
% x=m2(:,1);y=m2(:,2);z=m2(:,3);
% on=[1516 3755 6281 8986 11270]; % wire iron d
% off=[370 2917 5298 8061 10700]; %
% on=[2340 5131 7630 10270 12780]; % wire 100 20 60
% off=[1800 4560 6960 9536 12180]; %
% on=[1569 4554 7132 11370 13800]; % wire 100 20 60 10a
% off=[730 3695 6492 10050 13140]; %
% on=[1300 4700 9100 11750]; % wire 100 70 100 10a
% off=[840 4176 8506 11240]; %
% on=2750;% 7 8
% off=1000;
% on=[1870 7719 12080 15670 19320 23350 27190 34510 39750 44310 48550 53410 57840 63160 68200 72420 76670 82150 87790 92000 95570 100700];% wire 150 230 20 10a
% off=[1430 6603 11210 14980 18870 22820 26710 33830 38920 43620 47570 52390 56820 62130 66890 71070 75500 80910 86090 90680 94970 10200];
% on=[3000 7062 10820 15430 19570 23880 28320 32500 36700 41210 44770 48980 53000 57610 62460 66000 70000 73610 77680 81610 89000 93310];% wire2 150 230 20 10a
% off=[1943 5967 9599 13950 8090 22430 26830 30800 35600 39730 43650 47618 51450 56130 61190 64950 68800 72180 76400 80370 88000 92000];
% on=[1800 5600 9100 13000 16600 19700 23500 26600 30000 34270 38000 41340 44350 47840 51550 55570 59180 62450 66100 70880 73850 77700];% wire3 150 230 20 10a3
% off=[1250 5088 8600 12500 16000 19200 22500 26000 29500 33800 37340 40800 43900 47230 50850 55000 58530 61900 65420 69750 73060 76840];
% on=[5800 8052 10260 12510 14460];% wire3 150 230 190 19a
% off=[5620 7900 10040 12310 14340];
% on=[873 3516 6540 9534 12540 15740 18510 21360 24380 37800 31270 35000];% wire3 150 iron 230 120 10a
% off=[622 3137 6154 9016 12000 15090 18100 21000 23900 27350 30750 34450];
% [~,k]=size(on);
% wd(k,:)=0;
% for i=1:k
%     wd(i)=sqrt((x(on(i))-x(off(i)))^2+(y(on(i))-y(off(i)))^2+(z(on(i))-z(off(i)))^2);
% end
% wd
% loglog(d,wd)
% axis equal
% hold on
% plot(d,wd)


% a1=mm42(6600)-mm42(6800) %wire d test car
% a2=mm42(7800)-mm42(8100)
% a3=mm42(9100)-mm42(9400)
% a4=mm42(10900)-mm42(11100);
% a5=mm42(11700)-mm42(12000);
% a6=mm42(13000)-mm42(13300);
% a7=mm42(14000)-mm42(14100);
% a8=mm42(15200)-mm42(15400);
% a9=mm42(16000)-mm42(16300);
% a10=mm42(17400)-mm42(17600);
% mc=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10]'
% deg=[112.5 168.75 236.25 281.25 326.25 371.25 427.5 483.75 528.75 585]';
% d0=2*pi*4.52607796436989;
% d=deg./360*d0
% loglog(d,mc)
% plot(d,mc)
% axis equal tight


% plot(smooth(mm4,20))%%6 25 chair 指向磁場
% plot(smooth(mm42,20))
% h=-90.4;b=[h h h h h h h h h h h h h]';
% a=[120  100 80 60 40 20 0 -20 -40 -60 -80 -100 -120]';
% quiver(a,b,x,z,5)


% M0=[2 0.3 0.09 0.04 0.01 0.008]*100;
% D0=[50 100 150 200 250 300];
% DD0=(50:300);
% MM0=interp1(D0,M0,DD0,'cubic');
% 
% M=[4364.4 3126.8 645.1669 162.3161 97.6075 51.5263 30.8555 20.9089 15.4960 10.9715 7.4834 6.1782 4.8156 3.324 2.1273 1.5309 1.2144 0.9057 0.8152];
% D=[13     15     30       50       60      75      90      100     115     130     150    160    175    200   230    250    275    300    315];
% DD=(13:.1:315);
% 
% MM=interp1(D,M,DD,'cubic');%%curve
% 
% 
% % M2=[4790.3 3744.6 1583.7 955.7854 729.8330 177.1261 105.3062 60.1686 55.4310 33.0299 24.6006 16.6474 11.7061 7.927];
% M2=[4969.4 3600.8 1549.1 930.3427 708.4758 173.5085 103.3218 58.7619 54.1431 32.3240 23.9856 16.2510 11.4397 7.7166];
% D2=[13     15     22     27       30       50       60       73      75      90      100     115     130     150];
% DD2=(13:.1:150);
% 
% MM2=interp1(D2,M2,DD2,'cubic');
% 
% 
% bM=[967.8697 178.6601 60.9344];
% bD=[27       50       73];
% bDD=(27:73);
% bMM=interp1(bD,bM,bDD,'cubic');
% 
% % hold on
% % plot(DD0,MM0,'r',DD,MM,'b',DD2,MM2,'g',bDD,bMM,'m2')
% % 
% % plot(D0,M0,'ro',D,M,'bo',D2,M2,'go',bD,bM,'mo2')
% % ;
% % plot(22,744.6687,'ko') 
% % set(gca, 'xtick', D2);
% % grid on	
% % legend('模擬曲線','實測曲線','實測曲線2','n block data');
% % hold off
% 
% 
% %%%%reports
% % loglog(DD,MM,D,M,'o')
% % % loglog(DD2,MM2,'g',D2,M2,'o')
% % D2=[13     15     22     27       30       50       60       73      75      90      100     115     130     150    160    175    200   230    250    275    300    315];
% % set(gca, 'xtick', D2);
% % M=sort(M);
% % set(gca, 'ytick', M);
% % grid on	
% % xlabel('距離 (cm)')
% % ylabel('磁場強度 (uT)')
% 
% 
% % M=[4903.1 744.66 144.48]; %wall
% % M2=[4724.9 1549.1 624.11];
% % D=[13     22     32.1];
% % DD=(13:.1:32.1);
% % MM=interp1(D,M,DD,'cubic');
% % MM2=interp1(D,M2,DD,'cubic');
% % hold on
% % set(gca,'XScale' ,'log' ,'YScale' ,'log' );
% % loglog(DD,MM,D,M,'o')  
% % loglog(DD,MM2,D,M2,'o')
% % set(gca, 'xtick', D);
% % M=[144.48 624.11 744.66 1549.1 4903.1];
% % set(gca, 'ytick', M);
% % grid on
% % xlabel('距離 (cm)')
% % ylabel('磁場強度 (uT)')
% % hold off
% 
% % M=[967.8697 178.6601 60.9344]; %papers
% % D=[27       50       73];
% % DD=(27:.1:73);
% % MM=interp1(D,M,DD,'cubic');
% % 
% % plot(DD,MM,D,M,'o')
% % D2=[27 50 73];
% % set(gca, 'xtick', D2);
% % grid on	
% % xlabel('距離 (cm)')
% % ylabel('磁場強度 (uT)')
% 
% %%%%%
% 
% % 27cm n       m=967.8697
% % 27cm paper   m=967.6859
% % 50cm n       m=178.6601
% % 50cm paper   m=178.4894
% % 73cm n       m=60.9344
% % 73cm paper   m=60.9172
% % 22cm wall  m=744.6687
% % 22cm wall2 m=743.9897
% nblock=[967.8697 178.6601 60.9344]';
% block =[967.6859 178.4894 60.9172]';
% % wall=[744.6687 743.9897]';
% Dblock=[27 50 73]';
% Dblockwall=[22];
% 
% % ftest13cm wall1 m=4903.1  %27.6v 5.08a
% % ftest13cm wall2 m=4906.1
% % ftest32.1cm wall1 m=144.4837
% % ftest32.1cm wall2 m=144.4554
% % ftest13cm nwall1 m=4724.9
% % ftest13cm nwall2 m=4738.7
% % ftest32.1cm nwall1 m=624.1121
% % ftest32.1cm nwall2 m=623.0401
% wall  =[4903.1 144.4837]';
% wall2 =[4906.1 144.4554]';
% nwall =[4724.9 624.1121]';
% nwall2=[4738.7 623.0401]';
% Dftest=[13 32.1]';
% 
