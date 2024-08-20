clear all

%X=load('tests bal2.txt');
% X=load('fit data\12 4 wood with mag2.txt');
XL=load('fit data\12 6 wood 6mags.txt');
X=XL(:,[11 12 13]);
TWuT=45; % (高斯G)需設定 % 台灣磁場強度背景值45uT (0.45高斯) 
equals = ''; % no constraints by default

[ center, radii, evecs, v, chi2 ] = ellipsoid_fit_new( X, equals );

[a1 a2]=size(X);
corr_center=zeros(a1,3); corr_center(:,1)=center(1,1); corr_center(:,2)=center(2,1); corr_center(:,3)=center(3,1);
corr_radii=zeros(a1,3); corr_radii(:,1)=radii(1,1); corr_radii(:,2)=radii(2,1); corr_radii(:,3)=radii(3,1);
Xnew=(X-corr_center)./corr_radii*TWuT;
% +corr_center;

corr_center_X=center(1,1); % 取得 
corr_center_Y=center(2,1); % 取得
corr_center_Z=center(3,1); % 取得
corr_radii_X=radii(1,1); % 取得
corr_radii_Y=radii(2,1); % 取得
corr_radii_Z=radii(3,1); % 取得

hold on
grid on
grid minor
axis equal tight
xlabel('X');
ylabel('Y');
zlabel('Z');
title('校正前後比對')

scatter3(X(:,1),X(:,2),X(:,3),'r.')
scatter3(Xnew(:,1),Xnew(:,2),Xnew(:,3),'b.')
scatter3(center(1),center(2),center(3),'ro')
grid on
grid minor
axis equal tight
% %%%%
% %close all

% XL=load('fit data\12 4 wood with mag2.txt');
X=XL(:,[14 15 16]);
TWuT=45; % (高斯G)需設定 % 台灣磁場強度背景值45uT (0.45高斯) 
equals = ''; % no constraints by default

[ center, radii, evecs, v, chi2 ] = ellipsoid_fit_new( X, equals );

[a1 a2]=size(X);
corr_center2=zeros(a1,3); corr_center2(:,1)=center(1,1); corr_center2(:,2)=center(2,1); corr_center2(:,3)=center(3,1);
corr_radii2=zeros(a1,3); corr_radii2(:,1)=radii(1,1); corr_radii2(:,2)=radii(2,1); corr_radii2(:,3)=radii(3,1);
Xnew2=(X-corr_center2)./corr_radii2*TWuT;

corr_center_X2=center(1,1); % 取得
corr_center_Y2=center(2,1); % 取得 
corr_center_Z2=center(3,1); % 取得
corr_radii_X2=radii(1,1); % 取得
corr_radii_Y2=radii(2,1); % 取得
corr_radii_Z2=radii(3,1); % 取得

%%%%
% 
scatter3(X(:,1),X(:,2),X(:,3),'g.')
scatter3(Xnew2(:,1),Xnew2(:,2),Xnew2(:,3),'y.')
scatter3(corr_center_X2,corr_center_Y2,corr_center_Z2,'go')
scatter3(0,0,0,'black')
legend('m1校正前','m1校正後','m1校正前圓心','m2校正前','m2校正後','m2校正前圓心','0,0,0');
hold off

%% 二號磁力儀
close all
mo=load('elecmags run\12 4\LOG21.TXT');
% mo=load('hweyi run 6m 12v.txt');
% m=xlsread('pass fly\1 data fitted.xlsx','LOG5','H1:J3969');
m(:,1)=mo(:,14);m(:,2)=mo(:,15);m(:,3)=mo(:,16);

% figure(1)
% hold on
% scatter3(m(:,1),m(:,2),m(:,3))  %M1 模擬圖
% axis equal tight

m(:,1) = (m(:,1) - corr_center_X2 )./ corr_radii_X2 *TWuT;
m(:,2) = (m(:,2) - corr_center_Y2 )./ corr_radii_Y2 *TWuT;
m(:,3) = (m(:,3) - corr_center_Z2 )./ corr_radii_Z2 *TWuT;
m=m*0.01;
% scatter3(m(:,1),m(:,2),m(:,3))  %M1 模擬圖
% hold off
% 
% n=460;
% if(n>0)
% m(1:n,:)=[];
% end

%% 一號磁力儀
mo2=load('elecmags run\12 4\LOG21.TXT'); 
% mo2=mo;
% m2=xlsread('pass fly\2 data fitted.xlsx','LOG5','H1:J3969');
m2(:,1)=mo2(:,11);m2(:,2)=mo2(:,12);m2(:,3)=mo2(:,13);

% figure(2)
% hold on
% scatter3(m2(:,1),m2(:,2),m2(:,3))  %M2 模擬圖
% axis equal tight

m2(:,1) = (m2(:,1) - corr_center_X )./ corr_radii_X *TWuT;
m2(:,2) = (m2(:,2) - corr_center_Y )./ corr_radii_Y *TWuT;
m2(:,3) = (m2(:,3) - corr_center_Z )./ corr_radii_Z *TWuT;
m2=m2*0.01;
% scatter3(m2(:,1),m2(:,2),m2(:,3))  %M2 模擬圖
% hold off

%%

hold on
mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);
mm42=sqrt(m2(:,1).^2+m2(:,2).^2+m2(:,3).^2);
% plot(m)
% plot(m2)
plot(mm42)
plot(mm4)
% plot(m(:,1))
% plot(m(:,2))
% plot(m(:,3))
% plot(m2(:,1))
% plot(m2(:,2))
% plot(m2(:,3))
% legend('一號','二號')
% hold off


%% %%%%%% elecmag %%%%%%
% figure(2)
% tt=mo2(:,38); %frequency
%  plot(tt)
t=mo2(:,1);

% d=[80 160 240 320 401 481 561];% 3 elec d
% d=[80 160 240 320 401 481 561 642 719 820 900];% 3 elec d 9f
d=[2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10]*100;% 10/2 5 elc d test 10m
smm=25;%移動平均數
% x=smooth(m2(:,1),smm);y=smooth(m2(:,2),smm);z=smooth(m2(:,3),smm);
x=smooth(m2(:,1),smm);y=smooth(m2(:,2),smm);z=smooth(m2(:,3),smm);
xyz=sqrt(x.^2+y.^2+z.^2);

%%% use hz
% off=[17 57 90 130 163 200 240];% 2f 15s %%%3 em
% off=[15 56 100 132 170 215 268];% 4f 42s
% off=[6 39 70 99 166 201 244];% 5f 20s
% off=[12 46 83 121 170 225 271];% 6f 20s
% off=[12 51 96 142 190 245 290];% 8f 17s
% off=[6 43 82 116 163 198 237];% 9f 15s
% off=[6 43 82 116 163 198 237 294 370 424 460];% 9f

% off=[12 52 87 147 187 232 277];% 0f 12s %%%1 em?
% off=[12 64 100 148 185 229 288];% 2f 18s
% off=[13 51 82 120 163 204 241];% 4f 8s 235hz
% off=[18 61 95];% 5f 10N 8s
% off=[6 60 111 170 270 325];% 6f 10N 60v 29s 1000
% off=[12 60 120 167 220 276];% 7f 10N 24v re 43.5s
% off=(off)*200;
% off=off+800;
% on=off-1500; %

% off=[111966 135805 150740 168599 183109 196570 207754 218110 230072 243664 253950 269110 283709];
% off=off+800; % 8 27 hole test
% on=off-2000; 

% on=195000; %%9 25
% off=196000;
% on=557100;
% off=556900;


% on=[21 81 143 179 216 251 283 310 340 371 405 435 466 498 530 560 593];   %%%%%%%%%%%%%%[0m] 32s 5.5 %%% 10/2 em d test 10m
% on=[19 50 86 116 146 180 225 260 289 323 352 384 423 461 523 562 604 640 677 718 762]; %[0.5m] 30s 5.5
% on=[24 61 103 136 169 206 244 275 310 340 370 401 456 484 514 543 577 607 634 682 716];%[1m] 110s 5.5
% on=[20 56 93 123 154 182 222 255 284 315 347 381 410 444 481 516 549 591 624 657 690]; %[1.5m] 50s
% on=[16 49 85 118 150 180 210 241 274 308 340 373 406 470 522 561 595 646 702 734 768]; %[2m] 142s %%% 10/3 5 em d test 10m
% on=[20 63 93 125 156 187 215 252 279 310 347 383 413 443 480 522 550 580 629 672 710]; %[3m] 143s
% on=[19 51 80 112 142 185 221 260 291 326 357 390 419 493 526 572 617 665 739 800 843]; %[4m] 51s
% on=[15 48 83 118 148 180 213 249 282 319 350 401 436 475 519 564 606 649 682 728 795]; %[5m] 92s
% on=[23 60 95 134 167 195 231 265 316 454 487 521 561 594 634 679 717 748 791 846 889]; %[6m] 579s 

% on=[16 47 82 115 144 182 208 248 279 309 336 366 396 434 463 496 525 560 605 641 671]; %[7m] 110s 5.5 %%% 10/4 5 em d test 2,10m
% on=[15 46 80 111 154 187 201 251 285 326 356 396 424 453 486 543 575 613 660 704 739]; %[8m NG] s 
% on=[16 45 79 107 146 181 213 245 277 311 349 395 430 469 503 547 585 617 651 729 765]; %[8m] 580s 5.5 
% on=[15 60 97 140 193 227 262 293 347 408 445 484 518 572 619 662 696 730 790 816 862]; %[9m] 83s 5.5
% on=[15 46 77 113 147 178 208 250 316 382 428 468 504 539 575 615 654 685 746 797 838]; %[10m] 48s 5.5 
% on=[15 89 122 171 206 243 295 334 372 404 443 472 499 607 634 711 741 769 797 839]; %[20m] 114s 5
% on=(on+48)*200; 
% off=on-200*5.5;


% ticks=[on off];
% ticks=sort(ticks,'ascend');

%%% use time
%%% 1  2  3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21
% on=[0 200 305 405 495 590 680 780 875 970 1055 1165 1290 1385 1510 1665 2120 1815 1980 2470 2545 2620 2710 2780 2840 2900 3001 3075 3145 3210 3280 3350 3430 3520 3585 3650 3716.5 3790]+228; %%1~2 228s 5%%%10 9 glass 
% on=[0 70 125 180 250 305 360 410 475 525 580 635 695 750 805 859 925 995 1050 1290 1345 1400 1455 1510 1560 1610 1670 1720 1775 1825 1885 1935 1985 2040 2095 2150 2210 2260]+57;%%3~4 57s 5%%%10 9 glass 
% off=on-5;

% on=[0 102 143 188 231 270 305 363 411 456 495 537 613 660]+91.45;%10 17 oldwall
% off=on-5;

% on=[278 391 506 962 1294]+22;% 12 5m %11/14 tunnel d test
% on=[286 511 692 1199 1322 1621]+32;%13 7.4m(1x)-10m
% on=[149 245 810 893 1183 1302 1975 2078 2270 2346]+31;%14 10m(last2)-15m-20m 
% on=[241 323 567 653 1110 1215 1385 1459]+11;%15 17.5m-12.5m
% on=[71 157 210 271 336 392 451 518 583 643 705 766 825 886 957]+11;%11/19 star d test
% off=on+36;

% on=[32000 138500];
% off=[30000 140000];

% on=[51 7302.25]; %11/22 %%elecmag tm test
% off=on-13.5
% on=on+2.5*5;
% off=off+2.5*5;
% on=[0 2 7 13 21 29 41 56 79 121]*60+53.5; %11 25
% off=on+15;
% tem=[26.2 30.1 35.1 40.1 45.1 50.1 55.2 60 65 70];

% on=[35 85]-3.6; %11 26 LOG1
% off=[63 113]-3.6;

% on=[0 60 150 292 435 557 713 843 1028 1252 1403 1617 1795 2050 2323 2606 2973 3278 3780 4394 4910 5895 6674]+33; %11 26 LOG2
% off=on+30;

% on=[32 96]-3.6; %11 26 LOG1
% off=on+30;

% on=[93 167 285 363]; % 1 %11 28 
% on=[128 212 448 528]; % 2 %11 28 
% off=on+30;

% on=[33 114 222 317 404 490 573 725]+18;%11 29 %2 elecmag run
% off=on+45;

on=1;
off=2;

% ticks=[on off];
% ticks=sort(ticks,'ascend');
[~,k]=size(on);
ed(k,:)=0;
eds(k,3)=0;
ton(k,:)=0;
toff(k,:)=0;

for i=1:k
[~,ton(i)]=min(abs(t-on(i)));
[~,toff(i)]=min(abs(t-off(i)));
end
ticks=[ton toff];
% ticks(end+1)=250;
% ticks(end+1)=ton(1)+250+3350;
ticks=reshape (ticks, 1, numel(ticks));
ticks=sort(ticks,'ascend');


for i=1:k % use time
    ed(i)=sqrt((x(ton(i))-x(toff(i)))^2+(y(ton(i))-y(toff(i)))^2+(z(ton(i))-z(toff(i)))^2);
    eds(i,:)=[x(ton(i))-x(toff(i)) y(ton(i))-y(toff(i)) z(ton(i))-z(toff(i))]';
end
% for i=1:k % use hz
%     ed(i)=sqrt((x(on(i))-x(off(i)))^2+(y(on(i))-y(off(i)))^2+(z(on(i))-z(off(i)))^2);
%     eds(i,:)=[x(on(i))-x(off(i)) y(on(i))-y(off(i)) z(on(i))-z(off(i))]';
% end

% ed
% t(on(1))/60
% t(on(2))/60
% ed(i)=sqrt((xx1(i)-xx2(i))^2+(yy1(i)-yy2(i))^2+(zz1(i)-zz2(i))^2)
% eds
% loglog(ed)
% axis equal

% figure(1);
hold on
% % plot(d,ed)
% % plot(y)
% % plot(m2(:,1))
% % plot(m2(:,2))
% % plot(m2(:,3))
plot(x)
plot(y)
plot(z)
plot(mm42)
plot(xyz)

% G=0.001*(abs(smooth(mo2(:,8),100))+abs(smooth(mo2(:,9),100))+abs(smooth(mo2(:,10),100)));
% A1=0.01*(abs(smooth(mo2(:,2),100))+abs(smooth(mo2(:,3),100))+abs(smooth(mo2(:,4),100)));
% A2=0.01*(abs(smooth(mo2(:,5),100))+abs(smooth(mo2(:,6),100))+abs(smooth(mo2(:,7),100)));
G=[mo2(:,8) mo2(:,9) mo2(:,10)];
A1=[mo2(:,2) mo2(:,3) mo2(:,4)];
A2=[mo2(:,5) mo2(:,6) mo2(:,7)];
% plot(G)
% plot(A1)
% plot(A2)

set(gca, 'xtick', ticks);
grid on
xlabel('Time')
ylabel('Gauss')
hold off

%% wave magggg bro
% tags(k,2)=0;
% atag(k,:)=0;
% atag1(k,:)=0;
% waveN=10;
% magmin(waveN,k)=0;
% magmax(waveN,k)=0;
% 
% 
% hold on
% plot(xyz)
% filtCutOff = 0.9;
% [b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
% xyzL = filtfilt(b, a, xyz);
% plot(xyzL)
% for i=1:k
% testdata=xyz([ton(i):toff(i)],:);
% % testdata=xyz([toff(i):ton(i)],:);
% % tag=[250 200 270 160 140];%11/14
% % tag=[50 50 50 50 50 50];%
% % tag=[100 200 200 200 100 100 200 200 200 200];%
% % tag=[100 200 200 200 100 200 200 200];%
% 
% % tag=[100 300 150 300 150 450 350 350 350 350 450 450 550 150 150];%11/19
% 
% % tag=ones(size(on))*600;
% % tagg=ones(size(tag))*3500;
% % tag=ones(size(on))*50;
% % tagg=ones(size(tag))*3300;
% % tag=ones(size(on))*300; % 11 28 iron
% % tagg=ones(size(tag))*3300;
% tag=ones(size(on))*150; % 11 29 2elec mag
% tagg=ones(size(tag))*4400;
% 
% 
% % tagg=ones(size(tag))*3000; % 11 26 (on*1)
% 
% tags(:,1)=tag;
% tags(:,2)=tagg;
% 
% testdata(1:tags(i,1))=[];%
% testdata(tags(i,2):size(testdata))=[];%
% % figure(i);
% % hold on
% % plot(testdata)
% [b, a] = butter(1, (2*filtCutOff)/(1/0.01), 'low');
% testdata = filtfilt(b, a, testdata);
% % plot(testdata)
% c=findpeaks(testdata);
% IndMin=find(diff(sign(diff(testdata)))>0)+1;   %獲得局部最小值的位置
% IndMax=find(diff(sign(diff(testdata)))<0)+1;   %獲得局部最大值的位置
% % figure; hold on; box on;
% % plot(1:length(testdata),testdata);
% % plot(IndMin,testdata(IndMin),'r^')
% % plot(IndMax,testdata(IndMax),'k*')
% % legend('原數據','低通','波谷','波峰')
% 
% atag(i)=ton(i)+tags(i,1);
% atag1(i)=ton(i)+tags(i,1)+tags(i,2);
% % atag(i)=toff(i)+tags(i,1);
% % atag1(i)=toff(i)+tags(i,1)+tags(i,2);
% 
% IndMin(waveN+1:size(IndMin))=[]; % 使同長度
% IndMax(waveN+1:size(IndMax))=[];
% 
% magmin(:,i)=IndMin+atag(i);
% magmax(:,i)=IndMax+atag(i);
% plot(magmin(:,i),xyz(magmin(:,i)),'r^')
% plot(magmax(:,i),xyz(magmax(:,i)),'k*')
% end
%     
% plot(x)
% plot(y)
% plot(z)
% % plot(G)
% % plot(A1)
% % plot(A2)
% 
% tag=[atag atag1];
% tag=reshape (tag, 1, numel(tag));
% tag=sort(tag,'ascend');
% set(gca, 'xtick', tag);
% grid on
% hold off
% 
% set(gca, 'xtick', ticks);
% grid on
% xlabel('Time')
% ylabel('Gauss')
% % legend('原數據','低通','波谷','波峰','','','','','','','','','','','','','','','x','y','z')
% 
% ew(k,:)=0;
% ews(k,3)=0;
% wx(waveN,k)=0;
% wy(waveN,k)=0;
% wz(waveN,k)=0;
% for i=1:k % wave ,use time
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
% end
% ew
% ews
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
%% %%%%核康mag ? %%%%
% m(1:361,:)=[]; %361 12v %2031 3.3v
% m(1:8540,:)=[];m2(1:8540,:)=[];
% m(1:5421,:)=[];m2(1:5421,:)=[];


% t0=400;t1=1000; %背景
% t2=2000;t3=2700; %開磁鐵區間
% x0=sum(m([t0:t1],1))/(t1-t0); %取平均  D test
% y0=sum(m([t0:t1],2))/(t1-t0);
% z0=sum(m([t0:t1],3))/(t1-t0);
% x1=sum(m([t2:t3],1))/(t3-t2);
% y1=sum(m([t2:t3],2))/(t3-t2);
% z1=sum(m([t2:t3],3))/(t3-t2);
% 
% x20=sum(m2([t0:t1],1))/(t1-t0);%取平均  D test
% y20=sum(m2([t0:t1],2))/(t1-t0);
% z20=sum(m2([t0:t1],3))/(t1-t0);
% x21=sum(m2([t2:t3],1))/(t3-t2);
% y21=sum(m2([t2:t3],2))/(t3-t2);
% z21=sum(m2([t2:t3],3))/(t3-t2);%
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
% mm2=sqrt((m(:,1)-x0).^2+(m(:,2)-y0).^2+(m(:,3)-z0).^2);
% mm22=sqrt((m2(:,1)-x20).^2+(m2(:,2)-y20).^2+(m2(:,3)-z20).^2);
% mm4=sqrt(m(:,1).^2+m(:,2).^2+m(:,3).^2);
% mm42=sqrt(m2(:,1).^2+m2(:,2).^2+m2(:,3).^2);
% plot(m)
% plot(m2)
% plot(mm22)
% plot(mm2)
% plot(mm42)
% plot(mm4)
% plot(m(:,1))
% plot(m(:,2))
% plot(m(:,3))
% plot(m2(:,1))
% plot(m2(:,2))
% plot(m2(:,3))
% scatter3(m(:,1),m(:,2),m(:,3))
% scatter3(m2(:,1),m2(:,2),m2(:,3))
% legend('一號','二號')
% hold off
%%

% sm0=std(mm2([t0:t1],1));
% sm1=std(mm2([t2:t3],1));
% sm=[sm0 sm1]';
% rotor=(mo(:,17)+mo(:,18))/2;
% t=mo(:,1);
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
% ax=mo(:,2)*9.8;
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
% plot(mm42)
% plot(xyz)
% set(gca, 'xtick', ticks);
% grid on	
%% hall
% % H1=mo2(:,14)/1000;
% % H2=mo2(:,15)/1000;
% % H1=mo2(:,14);
% % H2=mo2(:,15);
% H1=smooth(mo2(:,14)-500,1);
% H2=smooth(mo2(:,15)-500,1);
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
% % plot(mm42)
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
% % plot(mm42)
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
%%%%%%%%%%% wire m d test%%%%%%%%%%
% d=[20 30 40 50 60 70 80 90 100 110 120 130 140 150]; %%  wire~~~~
% d=[230 220 210 200 190 180 170 160 150 140 130 120 110 100 90 80 70 60 50 40 30 20]; %%  wire~~~~
% x=m(:,1);y=m(:,2);z=m(:,3);
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


% a1=mm4(6600)-mm4(6800) %wire d test car
% a2=mm4(7800)-mm4(8100)
% a3=mm4(9100)-mm4(9400)
% a4=mm4(10900)-mm4(11100);
% a5=mm4(11700)-mm4(12000);
% a6=mm4(13000)-mm4(13300);
% a7=mm4(14000)-mm4(14100);
% a8=mm4(15200)-mm4(15400);
% a9=mm4(16000)-mm4(16300);
% a10=mm4(17400)-mm4(17600);
% mc=[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10]'
% deg=[112.5 168.75 236.25 281.25 326.25 371.25 427.5 483.75 528.75 585]';
% d0=2*pi*4.52607796436989;
% d=deg./360*d0
% loglog(d,mc)
% plot(d,mc)
% axis equal tight


% plot(smooth(mm42,20))%%6 25 chair 指向磁場
% plot(smooth(mm4,20))
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
% % plot(DD0,MM0,'r',DD,MM,'b',DD2,MM2,'g',bDD,bMM,'m')
% % 
% % plot(D0,M0,'ro',D,M,'bo',D2,M2,'go',bD,bM,'mo')
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
% % 27cm n       m2=967.8697
% % 27cm paper   m2=967.6859
% % 50cm n       m2=178.6601
% % 50cm paper   m2=178.4894
% % 73cm n       m2=60.9344
% % 73cm paper   m2=60.9172
% % 22cm wall  m2=744.6687
% % 22cm wall2 m2=743.9897
% nblock=[967.8697 178.6601 60.9344]';
% block =[967.6859 178.4894 60.9172]';
% % wall=[744.6687 743.9897]';
% Dblock=[27 50 73]';
% Dblockwall=[22];
% 
% % ftest13cm wall1 m2=4903.1  %27.6v 5.08a
% % ftest13cm wall2 m2=4906.1
% % ftest32.1cm wall1 m2=144.4837
% % ftest32.1cm wall2 m2=144.4554
% % ftest13cm nwall1 m2=4724.9
% % ftest13cm nwall2 m2=4738.7
% % ftest32.1cm nwall1 m2=624.1121
% % ftest32.1cm nwall2 m2=623.0401
% wall  =[4903.1 144.4837]';
% wall2 =[4906.1 144.4554]';
% nwall =[4724.9 624.1121]';
% nwall2=[4738.7 623.0401]';
% Dftest=[13 32.1]';
% 
