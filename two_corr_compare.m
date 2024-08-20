%% 橢球校正
clear all
close all

XL=load('fit data\4 15 proj new2b with hall without M.txt');  %corr1
TWuT=45; % (高斯G)需設定 % 台灣磁場強度背景值45uT (0.45高斯) 
equals = ''; % no constraints by default
X=XL(:,[11 12 13]); % mag1
[center, radii, evecs, v, chi2 ,corr_center ,corr_radii] = ellipsoid_fit_new( X, equals );
Xnew=(X-corr_center)./corr_radii*TWuT;
corr_center_X=center(1,1);
corr_center_Y=center(2,1);
corr_center_Z=center(3,1);
corr_radii_X=radii(1,1);
corr_radii_Y=radii(2,1);
corr_radii_Z=radii(3,1);


fitdata=xlsread('fit data\9 30 fit data.xlsx');   %corr2
Ba=fitdata(1:3,1);Ta=fitdata(1:3,2:4);Ka=fitdata(1:3,5:7);
Bg=fitdata(5:7,1);Tg=fitdata(5:7,2:4);Kg=fitdata(5:7,5:7);
Bm=fitdata(9:11,1);Tm2a=fitdata(9:11,2:4);


comparedata=load('C:\Users\henry chen\Desktop\magtest\proj\9 24 中油 iron pull big c\log1.csv');
% comparedata=load('proj\10 13 two corr compare\log15.csv'); %input
comparedata=comparedata(288000:298000,:);
datam=comparedata(:,11:13);
datamm4=sqrt(datam(:,1).^2+datam(:,2).^2+datam(:,3).^2);

m=comparedata(:,11:13);
m(:,1) = (m(:,1) - corr_center_X )./ corr_radii_X *TWuT; %corr1
m(:,2) = (m(:,2) - corr_center_Y )./ corr_radii_Y *TWuT;
m(:,3) = (m(:,3) - corr_center_Z )./ corr_radii_Z *TWuT;
mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);


data_corr=datam; %corr2
acc=comparedata(:,2:4); 
gyr=comparedata(:,8:10); 
acc2=acc;
gyr2=gyr;
[a1,~]=size(data_corr);
for i=1:a1   %校正數據
    acc2(i,:)=(Ta*Ka*(acc(i,:)'+Ba))';
    gyr2(i,:)=(Tg*Kg*(gyr(i,:)'+Bg))';
    data_corr(i,:)=(Tm2a*(comparedata(i,11:13)'+Bm))';
end
data_corr(:,1)=data_corr(:,1).*-1;
data_corr(:,2)=data_corr(:,2);
data_corr(:,3)=data_corr(:,3).*-1;
acc2(:,1)=acc2(:,1);
acc2(:,2)=acc2(:,2);
acc2(:,3)=acc2(:,3);
gyr2(:,1)=gyr2(:,1);
gyr2(:,2)=gyr2(:,2)*-1;
gyr2(:,3)=gyr2(:,3)*-1;

Xnew2=data_corr;
Xneww2=sqrt(Xnew2(:,1).^2+Xnew2(:,2).^2+Xnew2(:,3).^2);

corr1(1,1)=max(m(:,1));
corr1(2,1)=max(m(:,2));
corr1(3,1)=max(m(:,3));
corr1(1,2)=min(m(:,1));
corr1(2,2)=min(m(:,2));
corr1(3,2)=min(m(:,3))

corr2(1,1)=max(Xnew2(:,1));
corr2(2,1)=max(Xnew2(:,2));
corr2(3,1)=max(Xnew2(:,3));
corr2(1,2)=min(Xnew2(:,1));
corr2(2,2)=min(Xnew2(:,2));
corr2(3,2)=min(Xnew2(:,3))
%% m
close all

figure(1)
hold on
plot(datamm4)
plot(Xneww2)
plot(mm4)
hold off

figure(2)
hold on
plot(datam)
plot(Xnew2)
plot(m)
hold off
%% plot 三軸
close all
g=9.80665;j=1;k=1;w=1;
n=3;
for i=n:n+8
    figure(i)
    if(i<n+3)
        hold on
        plot(acc(:,j))
        plot(acc2(:,j)/g)
        hold off
        j=j+1;
    end
    if(i>n+2&&i<n+6)
        hold on
        plot(gyr(:,k))
        plot(gyr2(:,k)/pi*180)
        hold off
        k=k+1;
    end
    if(i>n+5)
        hold on
        plot(m(:,w))
        plot(Xnew2(:,w))
        hold off
        w=w+1;
    end
end
blue='old'
orange='new'
%% motion
hold on
scatter3(m(:,1),m(:,2),m(:,3))
scatter3(Xnew2(:,1),Xnew2(:,2),Xnew2(:,3))
hold off