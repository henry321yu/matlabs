clear; close all;


% fitdata=xlsread('9 30 fit data.xlsx');
% Ba=fitdata(1:3,1);Ta=fitdata(1:3,2:4);Ka=fitdata(1:3,5:7);
% Bg=fitdata(5:7,1);Tg=fitdata(5:7,2:4);Kg=fitdata(5:7,5:7);
% Bm=fitdata(9:11,1);Tm2a=fitdata(9:11,2:4);

%%% 橢球校正

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

%% old Corr_Mag
close all

mo1=load('C:\Users\henry chen\Desktop\magtest\proj\9 24 中油 iron pull big c\log0.csv');

mo=mo1;

[No,~]=size(mo); 
t0=[1 No];
smk=25; %移動平均係數
smk=1; %移動平均係數
testn='1'; %檔名
t=mo(:,1);

mo=mo(:,[1:4,8:13]);

m=mo(:,[8:10]);%% 一號磁力儀
m(:,1) = (m(:,1) - corr_center_X )./ corr_radii_X *TWuT;
m(:,2) = (m(:,2) - corr_center_Y )./ corr_radii_Y *TWuT;
m(:,3) = (m(:,3) - corr_center_Z )./ corr_radii_Z *TWuT;
m=m*0.01; %to Guess

%%% acc corr fit
corr_center_aX=0.0171511511529850;
corr_center_aY=-0.00318848470804418;
corr_center_aZ=-0.0142065034430159;
corr_radii_aX=1.03431459019208;
corr_radii_aY=1.01712801230964;
corr_radii_aZ=1.01255711628185;

acc1(:,1) = (mo(:,2) - corr_center_aX )./ corr_radii_aX;
acc1(:,2) = (mo(:,3) - corr_center_aY )./ corr_radii_aY;
acc1(:,3) = (mo(:,4) - corr_center_aZ )./ corr_radii_aZ;

mo(:,2:4) = acc1(:,1:3);
mo(:,8:10) = m(:,1:3);

ax=acc1(:,1);
ay=acc1(:,2);
az=acc1(:,3);
aGo=sqrt(mo(:,2).^2+mo(:,3).^2+mo(:,4).^2);
aG=sqrt(ax.^2+ay.^2+az.^2);
%acc1、acc2、mag1 Calibrated

GG=sqrt(mo(:,8).^2+mo(:,9).^2+mo(:,10).^2)*0.01;

mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);

x=smooth(m(:,1),smk);
y=smooth(m(:,2),smk);
z=smooth(m(:,3),smk);
smedm=[x y z];
xyz=sqrt(x.^2+y.^2+z.^2);
%% differend pull in big data(but separated)
mooton=[];
mootoff=[];
data=mo1;
mooton=[8.8,11.8,31.6,35.7,55,58.1,63.9,67,73.8,87.4,109.6]*10^4; %mo1 %9 24 中油 iron pull
mootoff=[10,13,33,36.6,56,59.15,64.8,67.8,74.6,88.5,110.8]*10^4;
% data=mo2;
% mooton=[2 4.4 8.6 11.3 14.3 17.2 19.9 22.4 25.4]*10^4; %mo2
% mootoff=[2.9 5.4 9.45 12.3 15.4 18.2 21 23.4 26.5]*10^4;
% data=mo3;
% mooton=[1,1.22,1.42,1.65,1.87,2.16,3.2,3.42,3.67,3.92,4.16,4.42,6.05,6.32,6.53,6.75]*10^5; %mo3
% mootoff=[1.09,1.3,1.52,1.73,1.97,2.26,3.3,3.5,3.76,4.01,4.26,4.53,6.13,6.41,6.63,6.85]*10^5;
% data=mo4;
% mooton=[5.22,5.44,6.36,6.58,6.77,7.06,7.5,8.13,10,10.57,10.9]*10^5; %mo4
% mootoff=[5.32,5.52,6.48,6.68,6.86,7.16,7.62,8.22,10.1,10.66,11]*10^5;
% data=mo5;
% mooton=[3.4 6.3 9.6 16 28 31.1 35.3 62 64.8 67.8]*10^4; %mo5
% mootoff=[4.2 7.1 10.4 17.2 28.8 32 36.2 62.9 65.8 68.8]*10^4;

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
fid=fopen(['debug\',char(filename)],'w');%寫入檔案路徑
pulldata=[t(pton(i):ptoff(i)) mo((pton(i):ptoff(i)),8:10) mo((pton(i):ptoff(i)),2:4) mo((pton(i):ptoff(i)),5:7)];
[r,c]=size(pulldata);            % 得到矩陣的行數和列數
 for w=1:r
  for v=1:c
  fprintf(fid,'%f\t',pulldata(w,v));
  end
  fprintf(fid,'\r\n');
 end
fclose(fid);
end
%% threshold filter waves (for pipe measuring)
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