function SF_Accel_Gyro_Live()
close all;
if (~isempty(instrfind))
    fclose(instrfind);
end
clear all;
s1 = serial('COM3', 'BaudRate', 115200);
time = datestr(now,0);
global f;
global time_step;
global sensor_model;
global Qd;
global P;
global sv;
global counter;
global accelUpdateFrequency;
global velocity_history;
global position_history;
global roll_cov_history;
global q_est_history;
global roll_error;
global gyro_bias_history;
global gyro_cov_history;
global accel_history;

roll_cov_history = [];
q_est_history = [];
roll_error = [];

f = figure;
sensor_model = [];
time_step = 1/50;
velocity_history = [];
position_history = [];
q_est_history = [];
counter = 0;
accelUpdateFrequency = 3;
accel_history = [];
accel_prev = [0 0 0]';
sensor_model = init_sensor_model(sensor_model, time_step);
sv.gyro_bias_est = [0 0 0];
sv.accel_bias_est = [0 0 0]';
sv.q_est = [];
sv.dcm_est = [];
sv.velocity_est = [0 0 0]';
sv.position_est = [0 0 0]';

Qd =  [sensor_model.Q_q    zeros(3,3);
    zeros(3,3)          sensor_model.Q_gyro_bias];

P = [sensor_model.P_q    zeros(3,3);
    zeros(3,3)          sensor_model.P_gyro_bias];

tic;
s1.BytesAvailableFcn = {@mycallback, @draw_rotating_cube}

fopen(s1);
%if (counter > 5000)
w = waitforbuttonpress;
%if w == 1
%s1.BytesAvailableFcn = {[]};
save('log.mat', 'position_history', 'velocity_history', 'q_est_history');
fclose(s1);
figure
plot(3*sqrt(roll_cov_history), 'r');
hold on;
plot(roll_error - 0.0396, 'b'); % 0.0.25 is the roll angle offset for my test set up
plot(-3*sqrt(roll_cov_history), 'r');
title({'Roll angle residual plot', 'Measurement Update Frequency = 1/3*Sampling Frequency'...
    'Qd = gyro_noise_covariance*sample time^2'}, 'interpreter', 'none', 'FontWeight','Normal');
ylabel('roll error sigma (radians)');
xlabel('time');
legend('3 sigma (roll)', 'roll angle residual', '-3 sigma (roll)');
ylim([-0.1, 0.1]);
var(roll_error)

figure;
plot(3*sqrt(gyro_cov_history(:,1)), 'r');
hold on;
plot(-3*sqrt(gyro_cov_history(:,1)), 'r');
plot((gyro_bias_history(:,1))-0.01, 'b');
legend('3 sigma (gyro bias(x))', 'gyro bias(x) error', '-3 sigma (gyro bias(x))');
ylabel('gyro bias(x) error');
xlabel('time');
%end

% plot position
figure
plot3(position_history(:,1), position_history(:,2), position_history(:,3));
pos = [0 0 0];
vel = [0 0 0];
for idx = 1: length(accel_history)
    pos(end+1,:) = pos(end,:) + vel(end,:)*time_step + 0.5*accel_history(idx,:)*time_step^2;
    vel(end+1,:) = vel(end,:) + accel_history(idx,:)*time_step;
end

% make running plot of orientation
vidObj = VideoWriter('rotating_cube.avi');
loops = length(q_est_history) - 1500;
F(loops) = struct('cdata',[],'colormap',[]);
open(vidObj);
for idx = 1500: length(q_est_history)
    draw_rotating_cube(position_history(idx, :), quat2dc(q_est_history(idx, :)), f);
  writeVideo(vidObj, getframe(gcf));
end
close(vidObj);
clear variables;
%end
end


function mycallback(obj, event, draw_rotating_cube)
% pause(1/250);
global f;
global time_step;
global sensor_model;
global Qd;
global P;
global sv;
global counter;
global accelUpdateFrequency;
global velocity_history;
global position_history;
global roll_cov_history;
global q_est_history;
global roll_error;
global gyro_bias_history;
global gyro_cov_history;
global accel_history;
%time_step = toc;
tic;
g = [0 0 1]';
GYRO_SCALE_FACTOR = 131;
ACCEL_SCALE_FACTOR = 16384;
DEGTORAD = pi/180;
raw_data = fscanf(obj, '%d %d %d %d %d %d %d');
if (length(raw_data) ~= 7) return; end
raw_data = raw_data';

scaled_data_gyro = raw_data(2:4)./GYRO_SCALE_FACTOR;
scaled_data_accel = raw_data(5:7)./ACCEL_SCALE_FACTOR;
gx = scaled_data_gyro(1)*DEGTORAD; % measured values
gy = scaled_data_gyro(2)*DEGTORAD;
gz = scaled_data_gyro(3)*DEGTORAD;

ax = scaled_data_accel(1);
ay = scaled_data_accel(2);
az = scaled_data_accel(3);

% disp(sprintf('%f %f %f\n', ax, ay, az));
pitch_est = asin(ax);
roll_est = -atan(ay/az);
% accel doesn't give heading estimate, assume 0
heading_est = 0;
% if this is the first time, initialize the quaternion
if (isempty(sv.q_est))
    sv.dcm_est = euler2dc(roll_est, pitch_est, heading_est)';
    sv.q_est = dc2quat(sv.dcm_est)';
    sv.gyro_bias_est = [scaled_data_gyro(1) scaled_data_gyro(2) scaled_data_gyro(3)]'*DEGTORAD;
end
% measured omega
omega_measured = [gx gy gz]';
accel_measured = [ax ay az]';

% estimated omega
omega_est = omega_measured - sv.gyro_bias_est;
accel_est = accel_measured - sv.accel_bias_est;

phi = omega_est*time_step;
vel = accel_est*time_step;
vel3 = make_skew_symmetric_3(vel);
phi3 = make_skew_symmetric_3(phi);

sv.q_est = apply_small_rotation(phi, sv.q_est);
sv.dcm_est = quat2dc(sv.q_est);
[r p y] = dc2euler(sv.dcm_est)
F = eye(6) +   [-phi3       eye(3)*time_step;
    zeros(3,6)];
%   Qd = diag([0.08^2 0.068^2 0.0482^2 0.0001^2 0.0001^2 0.0001^2])*time_step;
P = F*P*F' + Qd; % propagate covariance
% assuming constant acceleration model
v_orig = sv.velocity_est;
sv.velocity_est = sv.velocity_est + sv.dcm_est*vel - g*time_step;
v_final = sv.velocity_est;
sv.position_est = sv.position_est + ((v_orig + v_final)/2)*time_step;

% For logging purposes
velocity_history(end+1,:) = sv.velocity_est';
position_history(end+1,:) = sv.position_est';
gyro_bias_history(end+1,:) = sv.gyro_bias_est';
q_est_history(end+1, :) = sv.q_est';
accel_history(end+1, :) = sv.dcm_est*[ax ay az]';
counter = counter + 1;
if (mod(counter, accelUpdateFrequency) == 0)
    % apply accel measurements
    [sv, P] = process_accel_measurement_update2(sv, accel_est, P, sensor_model.R_accel, g);
end
roll_cov_history(end+1, :) = P(1,1);
gyro_cov_history(end+1, :) = [P(4,4) P(5,5) P(6,6)];
[r p y] = quat2eul(sv.q_est);
sv.gyro_bias_est;
roll_error(end+1, :) = r;
P_trace = trace(P(1:2,1:2))
end
