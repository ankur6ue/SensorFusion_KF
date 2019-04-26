function sim_gyro3()
% Adds image based measurements and position, velocity and accel bias to the state vector
close all;
time_step = 1/100;
image_update_frequency = 10;
global DEG2RAD;
DEG2RAD = pi/180;

% stores indicies and lengths of various elements of the state vector
global KF_SV_Offset;

KF_SV_Offset.pos_index 			= 1;
KF_SV_Offset.pos_length 		= 3;
KF_SV_Offset.vel_index 			= 4;
KF_SV_Offset.vel_length 		= 3;
KF_SV_Offset.orientation_index 	= 7;
KF_SV_Offset.orientation_length = 4;
KF_SV_Offset.gyro_bias_index 	= 11;
KF_SV_Offset.gyro_bias_length 	= 3;
KF_SV_Offset.accel_bias_index 	= 14;
KF_SV_Offset.accel_bias_length 	= 3;

global state_history;
state_history.len = 100;
state_history.yaw_est		= zeros(1,state_history.len);
state_history.pitch_est		= zeros(1,state_history.len);
state_history.roll_est		= zeros(1,state_history.len);
state_history.yaw_true		= zeros(1,state_history.len);
state_history.pitch_true	= zeros(1,state_history.len);
state_history.roll_true		= zeros(1,state_history.len);
state_history.omega_measured= zeros(1,state_history.len);
state_history.omega_true	= zeros(1,state_history.len);
state_history.state_uncertainity= zeros(1,state_history.len);
state_history.NESS			= zeros(1,state_history.len);
state_history.gyro_bias		= zeros(1,state_history.len);
state_history.accel_bias	= zeros(1,state_history.len);

num_sim_samples = 600;

initial_pos = [0 0 0];
initial_vel = [0 0 0];
initial_accel = [0 0 0];
g = [0 0 1]';

sim_state.num_sim_samples = num_sim_samples;
sim_state.trajectory_rotation = [0 0 0; 20 20 20; 20 0 0; 10 10 0;  30 20 20; 0 0 0]*DEG2RAD;
sim_state.trajectory_accel = [0 0 0; 0.25 0 0; -0.40 0.25 0; 0 -0.40 0.25; 0 0 -0.40; 0 0 0];
sim_state.velocity = zeros(sim_state.num_sim_samples, 3);
sim_state.position = zeros(sim_state.num_sim_samples, 3);
sim_state.accel = zeros(num_sim_samples, 3);
sim_state.orientation = zeros(num_sim_samples, 4);
sim_state.velocity(1,:) = initial_pos;
sim_state.position(1,:) = initial_vel;
sim_state.accel(1,:) = initial_accel;
sim_state.orientation(1,:) = eul2quat(sim_state.trajectory_rotation(1,1), sim_state.trajectory_rotation(1,2), sim_state.trajectory_rotation(1,3));

sim_state = init_simulator(sim_state, time_step, g);

global f4;

f1 = figure; % for running plot
f2 = figure; % for rotating cube
f3 = figure;
f4 = figure;

sensor_model = [];
sensor_model = init_sensor_model(sensor_model, time_step);

Qd = [sensor_model.Q_p  zeros(3,3)          zeros(3,3)          zeros(3,3)      zeros(3,3);
      zeros(3,3)        sensor_model.Q_v    zeros(3,3)          zeros(3,3)      zeros(3,3);
      zeros(3,3)        zeros(3,3)          sensor_model.Q_q    zeros(3,3)      zeros(3,3);
      zeros(3,3)        zeros(3,3)          zeros(3,3)          sensor_model.Q_gyro_bias zeros(3,3);
      zeros(3,3)        zeros(3,3)          zeros(3,3)          zeros(3,3)      sensor_model.Q_accel_bias];
  
P = [sensor_model.P_p   zeros(3,3)          zeros(3,3)          zeros(3,3)      zeros(3,3);
      zeros(3,3)        sensor_model.P_v    zeros(3,3)          zeros(3,3)      zeros(3,3);
      zeros(3,3)        zeros(3,3)          sensor_model.P_q    zeros(3,3)      zeros(3,3);
      zeros(3,3)        zeros(3,3)          zeros(3,3)          sensor_model.P_gyro_bias zeros(3,3);
      zeros(3,3)        zeros(3,3)          zeros(3,3)          zeros(3,3)      sensor_model.P_accel_bias];

ness = 0;

% Camera setup
% lets use a NED coordinate system.
% position of feature points in the global coordinate system
p_f_in_G = [10 -3 -2; 10 3 -2; 10 -3 2; 10 3 2];
J1_pfg = [];
numFeatures = length(p_f_in_G);
% std dev for the image pixel noise
sigma_image = 0.1;
imageW = 10000;
imageH = 10000;
fl_y = 1000;
fl_z = 1000;
K = [fl_y, 0, imageW/2; 0, fl_z, imageH/2]; % camera calibration matrix

image_measurements = [];

p_C_in_G = [0 0 0]';

% initial estimates
sv.position_est = sim_state.position(1,:)'; % col vectors
sv.velocity_est = sim_state.velocity(1,:)';
sv.gyro_bias_est = 0*[0.0127 0.0177 0.0067]';
sv.accel_bias_est = [0 0 0]';
q_prev = sim_state.orientation(1,:);
sv.q_est = q_prev';
sv.dcm_est = quat2dc(sv.q_est);
sv.p_C_in_G_est = sv.position_est;
sv.KF_SV_Offset = KF_SV_Offset;

% visualization variables
position_history = [];
position_true_history = [];
velocity_history = [];
q_est_history = [];
accel_bias_history = [];
gyro_bias_history = [];
P_history = [];
roll_error = [];
pitch_error = [];
yaw_error = [];

% For making a video
global vidObj1;
vidObj1 = VideoWriter('cube_trajectory.avi');
open(vidObj1);

for i = 2:num_sim_samples-1
	idx = i;
    q = sim_state.orientation(i,:);
    dcm_true = quat2dc(q);
    position_true = sim_state.position(i,:);
    p_C_in_G = position_true';
	% Obtain simulated IMU measurements
    [omega_measured, accel_measured, q_prev] = sim_imu_tick(sim_state, time_step, idx, sensor_model, q_prev);
    % Obtain simulated camera measurements. p_f_in_G contains the feature
    % locations and p_C_in_G the camera location in the global coordinate system
    % sigma_image is the std. dev of the pixel noise, K the camera
    % calibration matrix.
    image_measurements = sim_camera_tick(dcm_true, p_f_in_G, p_C_in_G, K, numFeatures, sigma_image);

    % Obtain estimated values by adjusting for current estimate of the
    % biases
    omega_est = omega_measured - sv.gyro_bias_est;
    accel_est = accel_measured - sv.accel_bias_est;
 
    % Multiply by time_delta to get change in orientation/position
    phi = omega_est*time_step;
    vel = accel_est*time_step;
    vel3 = make_skew_symmetric_3(vel);
    phi3 = make_skew_symmetric_3(phi);
    
	% Generate new estimates for q, position, velocity using equations of motion
	sv = update_state(sv, time_step, g, phi, vel);
    % State transition matrix
    F = eye(15) + [zeros(3,3) eye(3)*time_step zeros(3,3)   zeros(3,3)      zeros(3,3);
                   zeros(3,3) zeros(3,3)       sv.dcm_est*vel3 zeros(3,3)      -sv.dcm_est*time_step;
                   zeros(3,3) zeros(3,3)        -phi3       eye(3)*time_step zeros(3,3);
                   zeros(3,15);
                   zeros(3,15)];
 
    % Propagate covariance           
	P = F*P*F' + Qd; 
    position_history(end+1,:)       = sv.position_est';
    position_true_history(end+1,:)  = position_true';
    q_est_history(end+1,:)          = sv.q_est';
    accel_bias_history(end+1,:)     = sv.accel_bias_est';
    gyro_bias_history(end+1,:)      = sv.gyro_bias_est';
	P_history(end+1,:,:)            = P;
    
    if (mod(i, image_update_frequency) == 0)
        [sv, P] = process_image_measurement_update(sv, p_f_in_G, image_measurements, ...
            numFeatures, P, K, imageW, imageH, sigma_image); 
    end
    if (mod(i, 10000) == 0)
        % apply accel measurements
        [sv, P] = process_accel_measurement_update(sv, accel_est, P, sensor_model.R_accel, g);
	end
    P_trace = trace(P);
    [r1 p1 y1] = dc2euler(dcm_true); % true euler angles
    [r2 p2 y2] = dc2euler(sv.dcm_est); % noisy, biased angles
  %  show_running_plot(r1,p1,y1, r2, p2, y2, ness, omega_measured, omega, P_trace, f3);
    roll_error(end+1) = r2-r1; yaw_error(end+1) = y2-y1; pitch_error(end+1) = p2-p1;
   % pause(0.001);
    
   % draw_rotating_cube(sv.position_est, sv.dcm_est, f2);
   draw_rotating_cube(position_true', dcm_true, f2);
    pause(0.001);
end

close(vidObj1);

%figure;
%plot(position_history(:,1));
% solve for steady state P matrix
% H_accel = [make_skew_symmetric_3(g)*dcm_est zeros(3,3)];
% FF = F(7:12, 7:12);
% Qdd = Qd(7:12, 7:12);
% R = sensor_model.R_accel;
% [X,L,G] = dare(FF', H_accel', Qdd, R)
figure;
plot3(position_true_history(:,1), position_true_history(:,2), position_true_history(:,3));
hold on;
plot3(position_history(:,1), position_history(:,2), position_history(:,3));
xlabel('x (meters)'); ylabel('y (meters)'); zlabel('z (meters)');
legend('true position', 'estimated position');
P(7:9, 7:9)
P(1:3, 1:3)
figure;
plot(3*sqrt(P_history(:,1,1)), 'r')
hold on;
plot(-3*sqrt(P_history(:,1,1)), 'r')
plot(position_history(:,1) - position_true_history(:,1), 'b')
ylabel('Position(x) error sigma (m)');
xlabel('time');
legend('3 sigma Position(x) error (m)', 'position error residual', '-3 sigma Position(x) error(m)');
ylim([-0.5, 0.5]);

figure;
plot(3*sqrt(P_history(:,2,2)), 'r')
hold on;
plot(-3*sqrt(P_history(:,2,2)), 'r')
plot(position_history(:,2) - position_true_history(:,2), 'b')
ylabel('Position(y) error sigma (m)');
xlabel('time');
legend('3 sigma Position(y) error (m)', 'position error residual', '-3 sigma Position(y) error(m)');
ylim([-0.5, 0.5]);

figure;
plot(3*sqrt(P_history(:,7,7)), 'r')
hold on;
plot(-3*sqrt(P_history(:,7,7)), 'r')
plot(roll_error, 'b')
ylabel('roll error sigma (rad)');
xlabel('time');
legend('3 sigma roll error (rad)', 'roll error', '-3 sigma roll error (rad)');
ylim([-0.5, 0.5]);

sv.accel_bias_est
var(roll_error)
mean(roll_error)
mean(yaw_error)
mean(pitch_error)
%save('simulated_imu_data', 'gyro_measurements', 'accel_measurements');
end


