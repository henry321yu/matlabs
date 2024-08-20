%% 橢球校正
clear all

XL=load('fit data\4 15 proj new2b with hall without M.txt'); % fit data
mo=load('proj\9 24 中油 iron pull big c\log4.csv'); % measured data
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
% read
close all
% mo=load('C:\Users\henry chen\Desktop\LOG20.txt');
% mo=load('proj\7 2 中油 iron dic and pull\LOG55.txt');

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
xyz=sqrt(x.^2+y.^2+z.^2);
%% pull simulat
% close all
% msimulat=m(:,1)+0.044; % sample data
% logfeq=100;
% maxmag=0.5; %自設磁力量 gauss
% magfeq=2; % 頻率 hz
% ontime=0.04/0.46; % 開電佔比
% flag=14000; %star wait
% N=16000; %end
% magfeq=1/magfeq;
% for i=5:N
%     if(i>=flag)
%         msimulat(i:round(i+magfeq*logfeq*ontime),:)=msimulat(i:round(i+magfeq*logfeq*ontime),:)+maxmag;
%         msimulat(i-1,:)=msimulat(i-1,:)+maxmag*2/3;
%         msimulat(i-2,:)=msimulat(i-2,:)+maxmag/3;
%     flag=i+3+magfeq*logfeq;
%     end
% end
% hold on
% plot(msimulat)
% for i=5:No
%     if(abs(msimulat(i,1)-msimulat(i-2,1))>maxmag*0.9&&maxmag>0)
%         plot(i,msimulat(i,1),'r.')
%     end
% end
%% pull simulat2
close all
msimulat=m(:,1); % sample data
logfeq=400;
setmaxmag=0.0001; %自設磁力量 gauss
maxmag=[];
magfeq=2; % 頻率 hz
ontime=0.1/0.4; % 開電佔比
flag=620000; %star wait
flagi=flag;
tflag=7109; %自設最大磁力時間(s)
N=628000; %end
magfeq=1/magfeq;
tz=t(flag:N);
for i=5:N
    if(i>=flagi)
        maxmag(i)=(setmaxmag*abs(1-abs(tflag/abs(tz(i-flag+1)-tflag)))^0.33333);
        msimulat(i:round(i+magfeq*logfeq*ontime),:)=msimulat(i:round(i+magfeq*logfeq*ontime),:)+maxmag(i);
        msimulat(i-1,:)=msimulat(i-1,:)+maxmag(i)*2/3;
        msimulat(i-2,:)=msimulat(i-2,:)+maxmag(i)/3;
    flagi=i+3+magfeq*logfeq;
    i
    end
end
hold on
% plot(t(flag:N,:),msimulat(flag:N,:))
plot(maxmag)