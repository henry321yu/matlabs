function [w_complementary, phi_theta_psi_complementary]=Complementary_Filter(imu_data,a_parameter,dt_i,N,init_quat,InitializationTime)


%bias_w = mean(imu_data(1:InitializationTime,4:5),1);
%imu_data(:,4:5)=imu_data(:,4:5)-bias_w;
%bias_w = mean(imu_data(1:InitializationTime,4:6),1);
%imu_data(:,4:6)=imu_data(:,4:6)-bias_w;

%a_parameter=0.99;   % 越靠近1 越相信陀螺儀不會飄移
phi_accelerometer=zeros(N,1);
theta_accelerometer=zeros(N,1);
psi_accelerometer=zeros(N,1);
phi_complementary=zeros(N,1);
theta_complementary=zeros(N,1);
psi_complementary=zeros(N,1);
q=zeros(4,N);
quat_accelerometer=zeros(N,4);
w_temp=zeros(3,N);
for i=1:1                           % complementary filter 所累積的位態
    % Angle = a × (angle + gyro × dt) + (a -1) × (accelerometer)
    
    w=imu_data(i,4:6)/180*pi; % 徑度
    Ow = [0     -w(1)   -w(2)    -w(3);...
        w(1)   0       w(3)    -w(2);...
        w(2)  -w(3)    0        w(1);...
        w(3)   w(2)   -w(1)     0  ];
    
    q(:,i) = init_quat' + 0.5*Ow*dt_i(i,1)*init_quat' ;
    q(:,i) = q(:,i) / norm(q(:,i));
    [ Yaw2,Pitch2,Roll2 ]  = quat2angle(q(:,i)', 'ZYX');
    
    init_a = mean(imu_data(1:1,1:3),1);
    init_a = init_a / norm(init_a);
    
    psi_accelerometer_init =  Yaw2; % 徑度
 %   phi_accelerometer_init = atan2(init_a(2), init_a(3));  % 徑度
%    theta_accelerometer_init = -asin(init_a(1));  % 徑度
    
        theta_accelerometer_init =-atan2(init_a(1),sign(init_a(3))*(init_a(2)^2+init_a(3)^2)^0.5);  % 初始姿態 pitch
   phi_accelerometer_init =-atan2(-init_a(2),init_a(3));  % 初始姿態 roll
    
    
    quat_accelerometer_init = angle2quat(psi_accelerometer_init, theta_accelerometer_init, phi_accelerometer_init,  'ZYX' );
    
    q(:,i)= a_parameter * q(:,i) + (1 - a_parameter ) * (quat_accelerometer_init');
    
    q(:,i) = q(:,i) / norm(q(:,i));
    [ Yaw3,Pitch3,Roll3 ]  = quat2angle(q(:,i)', 'ZYX');
    phi_complementary(i,1)=Roll3;
    theta_complementary(i,1)=Pitch3;
    psi_complementary(i,1)=Yaw3;
    
    %        phi_complementary(i,1) .....
    %            = a_parameter * (init_phi + imu_data(i,4) * dt) + (1 - a_parameter ) * (phi_accelerometer(i,1)); % 度
    %        theta_complementary(i,1) .....
    %            = a_parameter * (init_theta + imu_data(i,5) * dt) + (1 - a_parameter) * (theta_accelerometer(i,1)); % 度
end

test_Roll2Pitch2Yaw2=zeros(N,3);
for i=2:N
    % Angle = a × (angle + gyro × dt) + (a -1) × (accelerometer)
    
    w=imu_data(i,4:6)/180*pi; % 徑度
    Ow = [0     -w(1)   -w(2)    -w(3);...
        w(1)   0       w(3)    -w(2);...
        w(2)  -w(3)    0        w(1);...
        w(3)   w(2)   -w(1)     0  ];
    
    q(:,i) = q(:,i-1) + 0.5*Ow*dt_i(i,1)*q(:,i-1) ;
    q(:,i) = q(:,i) / norm(q(:,i));
    [ Yaw2,Pitch2,Roll2 ]  = quat2angle(q(:,i)', 'ZYX');
    
    test_Roll2Pitch2Yaw2(i,:)=[Roll2,Pitch2,Yaw2]/pi*180;
    
    acc = imu_data(i,1:3);
    acc = acc / norm(acc);
    
 %   phi_accelerometer(i,1) = atan2(acc(2), acc(3));  % 徑度
 %   theta_accelerometer(i,1) = -asin(acc(1)); % 徑度
    
        theta_accelerometer(i,1) =-atan2(acc(1),sign(acc(3))*(acc(2)^2+acc(3)^2)^0.5);  % pitch
    phi_accelerometer(i,1) =-atan2(-acc(2),acc(3));  % roll
    
        psi_accelerometer(i,1) =  Yaw2; %phi_accelerometer(i,1)-(Roll2-Yaw2); % 徑度
        
    quat_accelerometer(i,:) = angle2quat(psi_accelerometer(i,1), theta_accelerometer(i,1), phi_accelerometer(i,1),  'ZYX' );
    
    if norm(quat_accelerometer(i,:)-q(:,i)')>norm(quat_accelerometer(i,:)+q(:,i)') && i>3
        quat_accelerometer(i,:)=-1*quat_accelerometer(i,:);  % 暫時刪掉 !!!!!!!!!!!!!!!
    end
    
    q(:,i)= a_parameter * q(:,i) + (1 - a_parameter ) * (quat_accelerometer(i,:)');
    

    a=[-q(2,i-1) -q(3,i-1) -q(4,i-1);
        q(1,i-1) -q(4,i-1) q(3,i-1);
        q(4,i-1) q(1,i-1) -q(2,i-1);
        -q(3,i-1) q(2,i-1) q(1,i-1)];

b=q(:,i)-q(:,i-1);
w_temp(:,i)=a\b*2/dt_i(i,1);




    q(:,i) = q(:,i) / norm(q(:,i));
    [ Yaw3,Pitch3,Roll3 ]  = quat2angle(q(:,i)', 'ZYX');
    phi_complementary(i,1)=Roll3; % 徑度
    theta_complementary(i,1)=Pitch3; % 徑度
    psi_complementary(i,1)=Yaw3; % 徑度
    
    %        phi_complementary(i,1) .....
    %        = a_parameter * (phi_complementary(i-1,1) + imu_data(i,4) * dt) + (1 - a_parameter ) * (phi_accelerometer(i,1)); % 度
    %        theta_complementary(i,1) .....
    %        = a_parameter * (theta_complementary(i-1,1) + imu_data(i,5) * dt) + (1 - a_parameter) * (theta_accelerometer(i,1)); % 度
end

phi_complementary=phi_complementary/pi*180;     % 度
theta_complementary=theta_complementary/pi*180;  % 度
psi_complementary=psi_complementary/pi*180;  % 度

phi_accelerometer=phi_accelerometer/pi*180;     % 度
theta_accelerometer=theta_accelerometer/pi*180;     % 度
psi_accelerometer=psi_accelerometer/pi*180;     % 度

phi_theta_psi_complementary=[phi_complementary theta_complementary psi_complementary];

phi_theta_psi_accelerometer=[phi_accelerometer theta_accelerometer psi_accelerometer];

complementary_rate=zeros(N,3);    % 欲求     % 由complementary filter 微分而得的角速度
for i=2:N
  %  quat=q(:,i);                  % 不確定是哪一種
    quat=(q(:,i)+q(:,i-1))/2;         % 不確定是哪一種
    q_star = [-quat(2) -quat(3) -quat(4);
        quat(1) -quat(4)  quat(3);
        quat(4)  quat(1) -quat(2);
        -quat(3)  quat(2)  quat(1)];
    complementary_rate(i,:)=q_star\(q(:,i)-q(:,i-1))/dt_i(i,1)/pi*180*2;
end
phi_complementary_rate=complementary_rate(:,1);
theta_complementary_rate=complementary_rate(:,2);
psi_complementary_rate=complementary_rate(:,3);

w_complementary=[phi_complementary_rate theta_complementary_rate psi_complementary_rate];

w_temp(1,:) = smooth(w_temp(1,:),5*60*200);
w_temp(2,:) = smooth(w_temp(2,:),5*60*200);
w_temp(3,:) = smooth(w_temp(3,:),5*60*200);

OVER=1;

end
