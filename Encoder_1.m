function [Encoder_X,turning_velocity,turning_degree] = Encoder_1( inputEncoder ,Encoder_threshold)
% clear all

% ��J�s�X���q�����GinputEncoder�A
% �Ѻ�X �t�׾���turning_velocity �P �ֿn�첾�q����Encoder_X

% inputEncoder = load('inputEncoder.txt'); % �s�X���q���ƾ� ���]�w
% mo=load('C:\Users\henry chen\Desktop\magtest\proj\2 6 hall\LOG20.txt');
% inputEncoder = [mo(:,1) mo(:,30)];
% Encoder_threshold=10; % �Ѽƶ��]�w
% circumference= 0.078671875; %���դϱ� % ���l�P��(����) ���]�w  %%%���c��
% circumference= 0.07853981633974483096156608458199; %�~�|��
circumference= 0.29; %���դϱ� % ���l�P��(����) ���]�w  %%% 3D�C�L�p��
% circumference= 0.305; %���դϱ� % ���l�P��(����) ���]�w  %%% �T�����

[N a2]=size(inputEncoder);

peak_max=inputEncoder(1,2);
peak_min=inputEncoder(1,2);
peak_max_i_temp=0;
peak_min_i_temp=0;

condition=0;
NO_peak_max=0;
NO_peak_min=0;
peak_max_i=zeros(N,1);
peak_min_i=zeros(N,1);
for i=1:N
    if condition==0 % ��}�l���ɭ�
        if peak_max<inputEncoder(i,2)
            peak_max=inputEncoder(i,2);
            peak_max_i_temp=i;
        end
        if peak_min>inputEncoder(i,2)
            peak_min=inputEncoder(i,2);
            peak_min_i_temp=i;
        end
        if peak_max-peak_min>Encoder_threshold......
                && inputEncoder(i,2)<inputEncoder(1,2)......
                && peak_max_i_temp>0 %2 21�[��
            condition=2;
            NO_peak_max=NO_peak_max+1;
            peak_max_i(NO_peak_max,1)=peak_max_i_temp;
            peak_max=peak_min;
        end
        if abs(inputEncoder(i,2)-peak_min)>Encoder_threshold......
                && inputEncoder(i,2)>inputEncoder(1,2)......
                && peak_min_i_temp>0 %2 21�[��
            condition=1;
            NO_peak_min=NO_peak_min+1;
            peak_min_i(NO_peak_min,1)=peak_min_i_temp;
            peak_min=peak_max;
        end
    end
    
    if condition==1 % �s�X�����ϤO�j�צb�����Ȥ��W ���b��M�i�p��
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
    
    if condition==2 % �s�X�����ϤO�j�צb�����Ȥ��U ���b��M�i����
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
for i=1:N % �ϤO�ѥ���t
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
for i=1:N % �ϤO�ѭt�ॿ
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
    if  turning_flag(i,1)==-0.5 || turning_flag(i+1,1)==-0.5 % �ϤO�ѭt�ॿ
        dis_degree(i,1)=turning_degree(i+1,1)-turning_degree(i,1); % ���׮t
    end
    if  turning_flag(i,1)==+0.5 || turning_flag(i+1,1)==+0.5 % �ϤO�ѥ���t
        dis_degree(i,1)=turning_degree(i,1)-turning_degree(i+1,1); % ���׮t
    end
end

for i=2:N-1
    ave_degree=(dis_degree(i-1,1)+dis_degree(i,1))/2; % �������׮t
    dis_time=(inputEncoder(i+1,1)-inputEncoder(i-1,1))/2; % �ɶ��t
    turning_velocity(i,1)=ave_degree/360*circumference/dis_time; % �t�׾���
    Encoder_X(i+1,1)=Encoder_X(i,1)+ave_degree/360*circumference; % �ֿn�첾�q����
end
% Encoder_X(i)

figure(1)
plot(inputEncoder(:,1),inputEncoder(:,2));
hold on
plot(inputEncoder(peak_max_i(:,1),1),inputEncoder(peak_max_i(:,1),2),'*');
hold on
plot(inputEncoder(peak_min_i(:,1),1),inputEncoder(peak_min_i(:,1),2),'+');
% figure(2)
% plot(inputEncoder(:,1),turning_degree(:,1));
% figure(3)
% plot(inputEncoder(:,1),turning_velocity(:,1));
% figure(4)
% plot(inputEncoder(:,1),Encoder_X(:,1));