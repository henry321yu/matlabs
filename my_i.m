%clear; close all; clc;
clear all
close all

global Tube_X0 Tube_Roll Tube_Pitch Mag Stick_XYZ;

InitializationTime=20000;  %6 ���Ϊ�l�R�m���ɶ��I�ƥ�
dt = 1/350; % 1/�[���W�v

%%%%% �ѼƳ]�w
IMU_ooriginal =load('C:\Users\henry chen\Desktop\temp\log347.csv');  %%%TEMP FOLDER%%%
% IMU_ooriginal = load('16475.csv'); % ���]�w ����
timeo=IMU_ooriginal(:,1);

% starT=900;  %�����ܫ��w�ɶ�   %50-90office 203-303corridor 210-303corridor 900-960corridor L
% [~,starTn]=min(abs(timeo-starT));
% starT2=960;
% [~,starTn2]=min(abs(timeo-starT2));
% IMU_original = IMU_ooriginal(starTn:starTn2,:);
IMU_original = IMU_ooriginal;

time=IMU_original(:,1);
[N,a1]=size(time); % N ������I��
start_n=1;



%uuu=1;

%if uuu==1

% ���{�����ت��� �Τ@�ڹq���K �Φh�ӺϤO�� �]�@�ӭy�� ��z�U�ӺϤO���Ҷq�����ϤO����

%IMU = load('1M 6mag 1�Ӻ��K �]�T�ӭy��_1.5��.txt'); % ���]�w
%    IMU_original = load('1 14 stair pvc�� �԰ʻP�w�I-3.txt'); % ���]�w ����

time=IMU_original(:,1);

[N,a1]=size(time); % N ������I��

imu_data=zeros(N,9); % IMU�[����� 1-3���[�t�� 4-6�����t�� 7-9���ϤO

imu_data(:,1)=IMU_original(:,2);  % 1~3���[�t��(�����V�U��1)
imu_data(:,2)=IMU_original(:,3);  % 1~3���[�t��
imu_data(:,3)=IMU_original(:,4); % 1~3���[�t��
imu_data(:,4)=IMU_original(:,5); % 4~6�����t�� ��쬰 ��/s  ???
imu_data(:,5)=IMU_original(:,6); % 4~6�����t�� ��쬰 ��/s
imu_data(:,6)=IMU_original(:,7); % 4~6�����t�� ��쬰 ��/s
imu_data(:,7)=IMU_original(:,8)/100;  % 7~9���ϤO(����)
imu_data(:,8)=IMU_original(:,9)/100;  % 7~9���ϤO(����)
imu_data(:,9)=-1*IMU_original(:,10)/100 *(-1);  % 7~9���ϤO(����)

acc=imu_data(:,1:3);
gyr=imu_data(:,4:6);

%% �����n������
gyr1=gyr;
gyr1=gyr1-mean(gyr1(1:InitializationTime,:));
sample=gyr1;
ix=cumtrapz(time,sample(:,1));
iy=cumtrapz(time,sample(:,2));
iz=cumtrapz(time,sample(:,3));
hold
plot(ix)
plot(iy)
plot(iz)
%%

% span=1;
% imu_data(:,1) = [smooth(imu_data(:,1),span)]'; % �`�N! �[�t�� ���Ƥ�  !!!!!!!!!!!!!
% imu_data(:,2) = [smooth(imu_data(:,2),span)]'; % �`�N! �[�t�� ���Ƥ�  !!!!!!!!!!!!!
% imu_data(:,3) = [smooth(imu_data(:,3),span)]'; % �`�N! �[�t�� ���Ƥ�  !!!!!!!!!!!!!

%%%% �p��ɶ��t
dt_i=zeros(N,1);
for i=2:N-1
    dt_i(i,1)=(time(i+1,1)-time(i-1,1))/2;
end
dt_i(1,1)=dt/2; dt_i(N,1)=dt/2;
%%%%%%%%%

for i=4:6 % �R�m�ɶq�쪺��Ƥ~�o�˥�
    %	imu_data(:,i)=4*sigma_gyro*randn(N,1)/180*pi;     % �Ȯɳo��
    imu_data(:,i)=imu_data(:,i)-mean(imu_data(1:InitializationTime,i)); % �R�m�ɶq�쪺��Ƥ~�o�˥�
    %  imu_data(:,i)=imu_data(:,i)*0; % �R�m�ɶq�쪺��Ƥ~�o�˥�
end


figure(121)
plot(time(:,1)-time(1,1),imu_data(:,1))% 1~3���[�t��
hold on
plot(time(:,1)-time(1,1),imu_data(:,2))
hold on
plot(time(:,1)-time(1,1),imu_data(:,3))
xlabel('�ɶ�(��)', 'FontSize',14); % X�y�жb
ylabel('�[�t��(g)', 'FontSize',14); % Y�y�жb
legend('X_b��V�[�t��','Y_b��V�[�t��','Z_b��V�[�t��');

figure(122)
plot(time(:,1)-time(1,1),imu_data(:,4))% 4~6�����t��
hold on
plot(time(:,1)-time(1,1),imu_data(:,5))
hold on
plot(time(:,1)-time(1,1),imu_data(:,6))
xlabel('�ɶ�(��)', 'FontSize',14); % X�y�жb
ylabel('���t��(��/��)', 'FontSize',14); % Y�y�жb
legend('X_b��V���t��','Y_b��V���t��','Z_b��V���t��');

figure(123)
plot(time(:,1)-time(1,1),imu_data(:,7))% 7~9���ϤO
hold on
plot(time(:,1)-time(1,1),imu_data(:,8))
hold on
plot(time(:,1)-time(1,1),imu_data(:,9))
xlabel('�ɶ�(��)', 'FontSize',14); % X�y�жb
ylabel('�ϤO�[����(����)', 'FontSize',14); % Y�y�жb
legend('M_X_b(X_b��V�ϤO���q)','M_Y_b(Y_b��V�ϤO���q)','M_Z_b(Z_b��V�ϤO���q)');

imu_data=imu_data';
init_a = mean(imu_data(1:3,1:InitializationTime),2);
init_a = init_a / norm(init_a);

init_psi =  0;   % ��l���A
init_theta =-atan2(init_a(1,1),sign(init_a(3,1))*(init_a(2,1)^2+init_a(3,1)^2)^0.5);  % ��l���A pitch
init_phi =-atan2(-init_a(2,1),init_a(3,1));  % ��l���A roll
init_quat = angle2quat(init_psi, init_theta, init_phi);
%%
close all;
a_parameter=0.99999;   % 0.999 �V�a��1 �V�۫H���������|�Ʋ�
[w_complementary, phi_theta_psi_complementary, phi_theta_psi_accelerometer]=Complementary_Filter(imu_data',a_parameter,dt_i,N,init_quat,InitializationTime);

figure(125)
plot(time(:,1)-time(1,1),phi_theta_psi_complementary(:,1))% ½�u���A
hold on
plot(time(:,1)-time(1,1),phi_theta_psi_complementary(:,2))% �ɨ����A
plot(time(:,1)-time(1,1),phi_theta_psi_complementary(:,3))% �ɨ����A
xlabel('�ɶ�(��)', 'FontSize',14); % X�y�жb
ylabel('���A(����)', 'FontSize',14); % Y�y�жb
legend('Roll(½�u��)','Pitch(�ɨ�)', 'FontSize',12);
%%
% quat = zeros(length(time), 4);
% AHRSalgorithm = AHRS('SamplePeriod',dt, 'Kp', 0, 'KpInit', 0);
% for i = 1:length(time)
%     AHRSalgorithm.UpdateIMU(deg2rad([gyrX(i) gyrY(i) gyrZ(i)]), [accX(i) accY(i) accZ(i)]);
%     quat(i,:) = AHRSalgorithm.Quaternion;
% end

%%
acc1=acc;
acc1G=acc1*0;
acc1G(:,3)=1;
accA=gyr;
fs = 1/dt; % Sampling rate

phi_theta_psi_complementary1=phi_theta_psi_complementary;
fc = 5; % Cut off frequency
[b,a] = butter(2,fc/(fs/2),'low'); % Butterworth filter of order 6
% phi_theta_psi_complementary1= filtfilt(b, a, phi_theta_psi_complementary);



sample=gyr;
ix=cumtrapz(time,sample(:,1));
iy=cumtrapz(time,sample(:,2));
iz=cumtrapz(time,sample(:,3));

for i=1:N
    accA(i,1)=atan2(acc1(i,2),sqrt(acc1(i,1)^2+acc1(i,3)^2));
    accA(i,2)=atan2(acc1(i,1),sqrt(acc1(i,2)^2+acc1(i,3)^2));
    accA(i,3)=atan2(sqrt(acc1(i,1)^2+acc1(i,2)^2),acc1(i,3));
end

Yaw=phi_theta_psi_complementary1(:,3)/180*pi;
Pitch=phi_theta_psi_complementary1(:,2)/180*pi;
Roll=phi_theta_psi_complementary1(:,1)/180*pi;
% Yaw=phi_theta_psi_accelerometer(:,3)/180*pi;
% Pitch=phi_theta_psi_accelerometer(:,2)/180*pi;
% Roll=phi_theta_psi_accelerometer(:,1)/180*pi;
% Yaw=w_complementary(:,3)/180*pi; ..??
% Pitch=w_complementary(:,2)/180*pi;
% Roll=w_complementary(:,1)/180*pi;

dcm_ZYX = angle2dcm(Yaw,Pitch,Roll,'ZYX' );
acc11=acc1;
gyr1=gyr;

% for i=1:N   %% accelerometer angle
%     pitA(i)=atan2(acc1(i,1),sqrt(acc1(i,2)^2+acc1(i,3)^2));
%     rolA(i)=atan2(acc1(i,2),sqrt(acc1(i,1)^2+acc1(i,3)^2));
%     yawA(i)=atan2(sqrt(acc1(i,1)^2+acc1(i,2)^2),acc1(i,3));
% end
%quat  %%%%

quat = angle2quat(Yaw,Pitch,Roll,'ZYX' );
dcm_ZYX =quat2dcm(quat);

% quata = angle2quat(yawA,pitA,rolA,'ZYX' );
% gyr1=deg2rad(gyr1);
% quatg = angle2quat(gyr1(:,3),gyr1(:,2),gyr1(:,1),'ZYX' );
%%%%
% temp=quatmultiply(quatg,quat);
% [gyr11(:,3) gyr11(:,2) gyr11(:,1)]=quat2angle(temp);
% gyr11=rad2deg(gyr11);

for i=1:N
    acc11(i,:)=dcm_ZYX(:,:,i)*acc11(i,:)';
end
for i=1:N
    acc1G(i,:)=dcm_ZYX(:,:,i)*acc1G(i,:)';
end

acc111=acc11-acc1G; %% acc ����O[0 0 1]

gyr1=gyr1-mean(gyr1(1:InitializationTime,:)); %���lbia

for i=1:N
    gyr1(i,:)=dcm_ZYX(:,:,i)*gyr1(i,:)';
end

close all;

figure(110)
subplot (2, 1,1)
plot(acc111)
hold
% plot(stationary)
subplot (2, 1,2)
plot(gyr1)
hold
% plot(stationary*10)

AG=sqrt(gyr(:,1).^2+gyr(:,2).^2+gyr(:,3).^2);

fc = 30; % Cut off frequency
[b,a] = butter(2,fc/(fs/2),'low');
AG = filtfilt(b, a, AG);
fc = 1; % Cut off frequency
[b,a] = butter(2,fc/(fs/2),'low');
AA1 = abs(filtfilt(b, a, acc111));
accA = abs(filtfilt(b, a, accA));
% fc = 0.001; % Cut off frequency
% [b,a] = butter(2,fc/(fs/2),'high');
% AA11 = abs(filtfilt(b, a, AA1));

AA=sqrt(AA1(:,1).^2+AA1(:,2).^2+AA1(:,3).^2);

stationary=AG<0.3;
stationaryA=AA<0.014;



fc = 1; % Cut off frequency
[b,a] = butter(2,fc/(fs/2),'low');
acc111= filtfilt(b, a, acc111);

% fc = 0.007; % Cut off frequency
fc = 0.002; % Cut off frequency
[b,a] = butter(2,fc/(fs/2),'high'); % Butterworth filter of order 6
% acc111 = filtfilt(b, a, acc111);

gyr11=gyr1;
% fc = 0.001; % Cut off frequency
fc = 0.001; % Cut off frequency
[b,a] = butter(2,fc/(fs/2),'high'); % Butterworth filter of order 6
% gyr11= filtfilt(b, a, gyr1);

% fc = 0.001; % Cut off frequency
fc = 1; % Cut off frequency
[b,a] = butter(2,fc/(fs/2),'low'); % Butterworth filter of order 6
% gyr11= filtfilt(b, a, gyr1);

figure(118)
subplot (2, 1,1)
plot(acc111)
hold
plot(stationaryA)
plot(AA*10)
subplot (2, 1,2)
plot(gyr11)
hold
plot(stationary*10)

c=0.98
CC= zeros(N,3);
CC(:,3)=gyr1(:,3);
CC(:,2)=gyr(:,2)*c + accA(:,2)*(1-c);
CC(:,1)=gyr(:,1)*c + accA(:,1)*(1-c);
%%
va= zeros(N,3);
vg= zeros(N,3);
vc= zeros(N,3);
for i = 2:N
    if(stationary(i) == 0)
        if(stationaryA(i)==0)
            va(i,:) = va(i-1,:) + acc111(i,:)*9.80665*(time(i)-time(i-1));
        end
        vg(i,:) = vg(i-1,:) + gyr11(i,:)*(time(i)-time(i-1));
        vc(i,:) = vc(i-1,:) + CC(i,:)*(time(i)-time(i-1));
    else
        va(i,:)=[0 0 0];
        vg(i,:) = vg(i-1,:) + gyr11(i,:)*(time(i)-time(i-1));
        vc(i,:) = vc(i-1,:) + CC(i,:)*(time(i)-time(i-1));
    end
end

% for i = 2:N
%         va(i,:) = va(i-1,:) + acc111(i,:)*9.80665*(time(i)-time(i-1));
%         vg(i,:) = vg(i-1,:) + gyr11(i,:)*(time(i)-time(i-1));
%     if(stationary(i) == 1)
%         va(i,:)=[0 0 0];
%     end
% 
% end


velDrift = zeros(size(va));
stationaryStart = find([0; diff(stationary)] == -1);
stationaryEnd = find([0; diff(stationary)] == 1);
for i = 1:numel(stationaryEnd)
    driftRate = va(stationaryEnd(i)-1, :) / (stationaryEnd(i) - stationaryStart(i));
    enum = 1:(stationaryEnd(i) - stationaryStart(i));
    drift = [enum'*driftRate(1) enum'*driftRate(2) enum'*driftRate(3)];
    velDrift(stationaryStart(i):stationaryEnd(i)-1, :) = drift;
end

% Remove integral drift �h���n���}��
v1=va;
v1=v1-velDrift;
figure(113)
plot(va)
figure(114)
plot(v1)


pa=zeros(N,3);
for i=2:N
    pa(i,:) = pa(i-1,:) + va(i,:)*(time(i)-time(i-1));    % integrate velocity to yield position �N�t�׾�X��}�A��m
end
p1=zeros(N,3);
for i=2:N
    p1(i,:) = p1(i-1,:) + v1(i,:)*(time(i)-time(i-1));    % integrate velocity to yield position �N�t�׾�X��}�A��m
end

figure(115)
plot(pa)
pa(9660,1)
p1(9660,1)
figure(116)
plot(p1)
close all;
figure(117)
subplot (2, 1,1)
plot(va)
subplot (2, 1,2)
plot(v1)

figure(118)
subplot (2, 1,1)
plot(pa(:,1))
hold
plot(p1(:,1))
subplot (2, 1,2)
plot(vg)
% hold
% plot(vc)

%%
pa(end,1)
pa(end,1)


gyr1=gyr;
gyr1=gyr1-mean(gyr1(1:InitializationTime,:));

figure(125)
hold on
plot(time(:,1)-time(1,1),phi_theta_psi_complementary(:,2))% ½�u���A
plot(time(:,1)-time(1,1),pitA)

cumtrapz(time,gyr1(:,2));
% HP filter accelerometer data HP�L�o���[�t�׭p�ƾ�
filtCutOff = 0.05;
[b, a] = butter(1, (2*filtCutOff)/(1/0.0025), 'high');
ans = filtfilt(b, a, ans);
cumtrapz(time,ans);
plot(time(:,1)-time(1,1),ans)



xlabel('�ɶ�(��)', 'FontSize',14); % X�y�жb
ylabel('���A(����)', 'FontSize',14); % Y�y�жb
legend('Com','acc', 'FontSize',12);
%%
% ~~~~~~~~~ ���ਤ�סB�|���� �P �ഫ�x�} �� �����ഫ �����
Yaw=yyaw(1000)/180*pi;
Pitch=ppit(1000)/180*pi;
Roll=rrol(1000)/180*pi;
dcm_ZYX = angle2dcm(  Yaw, Pitch, Roll,  'ZYX' );
% �NGLOBE���Шt���নLOCAL���Шt���L�{�̧Ǭ�~
% �V�W����Z , �V���f�ɰw���Y����Yaw
% �V��(�V����)����Y , �C�Y����Pitch
% �V�_(�V�e)����X , �Y�V�k�ɬ���Roll
% dcm_ZYX = angle2dcm(  Yaw, Pitch, Roll,  'ZYX' );
% �YGLOBE���Шt�Τ����@�V�q[Ga;Gb;Gc]
% �h�bLOCAL���Шt�Τ��ӦV�q��Ū�ƫh��[La;Lb;Lc]
% [La;Lb;Lc] = dcm_ZYX * [Ga;Gb;Gc]
%
% dcm_ZYX = angle2dcm(  Yaw, Pitch, Roll,  'ZYX' );
% [La;Lb;Lc] = dcm_ZYX * [Ga;Gb;Gc]
%
% quat = angle2quat( Yaw, Pitch, Roll,  'ZYX' );
% dcm = quat2dcm(quat);
% [La;Lb;Lc] = dcm * [Ga;Gb;Gc]

rx=[ 1  0   0;
    0   cos(Roll)  sin(Roll);
    0  -sin(Roll)  cos(Roll) ];
ry=[ cos(Pitch)  0  -sin(Pitch);
    0   1  0;
    sin(Pitch)  0  cos(Pitch); ];
rz=[ cos(Yaw)  sin(Yaw)  0;
    -sin(Yaw)   cos(Yaw)  0;
    0   0   1 ];

%  rx*ry*rz �����G���P��  angle2dcm(  Yaw, Pitch, Roll,  'ZYX' )
%  �Y�Y
dcm_ZYX = angle2dcm(  Yaw, Pitch, Roll,  'ZYX' );
dcm_ZYX2 = rx*ry*rz;
%  �h dcm_ZYX = dcm_ZYX2

dcm_ZYX = angle2dcm(  Yaw, Pitch, Roll,  'ZYX' );
[ Yaw2,Pitch2,Roll2 ] = dcm2angle( dcm_ZYX, 'ZYX' );
% �h[Yaw2,Pitch2,Roll2]�P[Yaw,Pitch,Roll]�۵�

quat = angle2quat( Yaw, Pitch, Roll,  'ZYX' );
[ Yaw2,Pitch2,Roll2 ] = quat2angle(quat, 'ZYX');
% �h[Yaw2,Pitch2,Roll2]�P[Yaw,Pitch,Roll]�۵�

dcm = quat2dcm(quat);
quat2 = dcm2quat(dcm);
% �h quat2 �P quat �۵�

% �Y     dcm_ZYX  = angle2dcm(      Yaw,       Pitch,       Roll,  'ZYX' );  % �Y�� C_n_b
% �B     dcm_ZYX2 = angle2dcm( -1*Roll,  -1*Pitch,  -1*Yaw,  'XYZ' );  % �Y�� C_b_n
% �h    dcm_ZYX2*dcm_ZYX �N���� I (���x�})
% P.S.  dcm_ZYX=dcm_ZYX2��  �B  dcm_ZYX=inv(dcm_ZYX2)
