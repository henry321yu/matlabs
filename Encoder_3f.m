function [Encoder_X,turning_velocity,turning_degree,averageV] = Encoder_3f( inputEncoder ,Encoder_threshold,circumference)
% clear all

% Load and preprocess input data
% mo = load('C:\Users\sgrc-325\Desktop\git\pipes\data\1101_1.CSV');  % Data file path
% moo = mo(389000:392000, :);  % Select specific range
% inputEncoder = [moo(:, 1), moo(:, 19)];

% Define parameters
% Encoder_threshold = 70; % Threshold for peak detection
% circumference = 0.29045;   % Wheel circumference in meters

% Initialize variables
[N, ~] = size(inputEncoder);
peak_max = inputEncoder(1, 2);
peak_min = inputEncoder(1, 2);
condition = 0;
NO_peak_max = 0;
NO_peak_min = 0;
peak_max_i = zeros(N, 1);
peak_min_i = zeros(N, 1);
% Detect peaks and valleys
for i = 1:N
    if condition == 0  % Initial state
        if peak_max < inputEncoder(i, 2)
            peak_max = inputEncoder(i, 2);
            peak_max_i_temp = i;
        end
        if peak_min > inputEncoder(i, 2)
            peak_min = inputEncoder(i, 2);
            peak_min_i_temp = i;
        end
        if peak_max - peak_min > Encoder_threshold && inputEncoder(i, 2) < inputEncoder(1, 2)
            condition = 2;
            NO_peak_max = NO_peak_max + 1;
            peak_max_i(NO_peak_max) = peak_max_i_temp;
            peak_max = peak_min;
        end
        if abs(inputEncoder(i, 2) - peak_min) > Encoder_threshold && inputEncoder(i, 2) > inputEncoder(1, 2)
            condition = 1;
            NO_peak_min = NO_peak_min + 1;
            peak_min_i(NO_peak_min) = peak_min_i_temp;
            peak_min = peak_max;
        end
    elseif condition == 1  % Searching for peaks
        if peak_max < inputEncoder(i, 2)
            peak_max = inputEncoder(i, 2);
            peak_max_i_temp = i;
        end
        if peak_min > inputEncoder(i, 2)
            peak_min = inputEncoder(i, 2);
            peak_min_i_temp = i;
        end
        if peak_max - peak_min > Encoder_threshold
            condition = 2;
            NO_peak_max = NO_peak_max + 1;
            peak_max_i(NO_peak_max) = peak_max_i_temp;
            peak_max = peak_min;
        end
    elseif condition == 2  % Searching for valleys
        if peak_max < inputEncoder(i, 2)
            peak_max = inputEncoder(i, 2);
            peak_max_i_temp = i;
        end
        if peak_min > inputEncoder(i, 2)
            peak_min = inputEncoder(i, 2);
            peak_min_i_temp = i;
        end
        if peak_max - peak_min > Encoder_threshold
            condition = 1;
            NO_peak_min = NO_peak_min + 1;
            peak_min_i(NO_peak_min) = peak_min_i_temp;
            peak_min = peak_max;
        end
    end
end

peak_max_i = peak_max_i(1:NO_peak_max);
peak_min_i = peak_min_i(1:NO_peak_min);

% Initialize arrays for calculations
turning_flag = zeros(N, 1);
turning_degree = zeros(N, 1);
turning_velocity = zeros(N, 1);
dis_degree = zeros(N, 1);
Encoder_X = zeros(N, 1);

turning_flag(peak_max_i) = 1;
turning_flag(peak_min_i) = -1;

% Calculate turning degree and velocity
turning_condition = 0;
for i = 1:N
    if turning_flag(i) == 1
        turning_condition = 1;
        temp1 = i;
    end
    if turning_condition == 1 && turning_flag(i) == -1
        temp2 = i;
        turning_degree(temp1:temp2) = (inputEncoder(temp1:temp2, 2) - inputEncoder(temp2, 2)) * 180 / (inputEncoder(temp1, 2) - inputEncoder(temp2, 2)) - 90;
        turning_flag(temp1+1:temp2-1) = 0.5;
    end
end

turning_condition = 0;
for i = 1:N
    if turning_flag(i) == -1
        turning_condition = 1;
        temp1 = i;
    end
    if turning_condition == 1 && turning_flag(i) == 1
        temp2 = i;
        turning_degree(temp1:temp2) = (inputEncoder(temp1:temp2, 2) - inputEncoder(temp1, 2)) * 180 / (inputEncoder(temp2, 2) - inputEncoder(temp1, 2)) - 90;
        turning_flag(temp1+1:temp2-1) = -0.5;
    end
end

% Calculate dis_degree and turning_velocity
for i = 1:N-1
    if turning_flag(i) == -0.5 || turning_flag(i+1) == -0.5
        dis_degree(i) = turning_degree(i+1) - turning_degree(i);
    end
    if turning_flag(i) == 0.5 || turning_flag(i+1) == 0.5
        dis_degree(i) = turning_degree(i) - turning_degree(i+1);
    end
end

for i = 2:N-1
    ave_degree = (dis_degree(i-1) + dis_degree(i)) / 2;
    dis_time = (inputEncoder(i+1, 1) - inputEncoder(i-1, 1)) / 2;
    turning_velocity(i) = ave_degree / 360 * circumference / dis_time;
%     Encoder_X(i+1) = Encoder_X(i) + ave_degree / 360 * circumference;
end

% Smooth turning_velocity
turning_velocity = smoothdata(turning_velocity, 'movmean', 100);

for i = 2:N-1
    dis_time = (inputEncoder(i+1, 1) - inputEncoder(i-1, 1)) / 2;
    Encoder_X(i+1) = Encoder_X(i) + turning_velocity(i)*dis_time;
end

averageV=[];
for i=1:N
    if(turning_velocity(i)>0)
        averageV(end+1)=turning_velocity(i);
    end
end

% Plot results
figure(1)
hold on
plot(inputEncoder(:, 1), inputEncoder(:, 2));
plot(inputEncoder(peak_max_i, 1), inputEncoder(peak_max_i, 2), '*');
plot(inputEncoder(peak_min_i, 1), inputEncoder(peak_min_i, 2), '+');
hold off
title('Hall Sensor');
% 
% figure(2)
% plot(inputEncoder(:, 1), turning_degree);
% title('Turning Degree');
% 
% figure(3)
% plot(inputEncoder(:, 1), turning_velocity);
% title('Turning Velocity');
% 
% figure(4)
% plot(inputEncoder(:, 1), Encoder_X);
% title('Position');
