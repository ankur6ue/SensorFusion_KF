function SF_Accel_Gyro()
% Adds image based measurements and position, velocity and accel bias to the state vector
%close all;
time_step = 1/50;

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

num_sim_samples = 1000;


initial_pos = [0 0 0];
initial_vel = [0 0 0];
initial_accel = [0 0 0];
g = [0 0 1]';

sim_state.num_sim_samples = num_sim_samples;
sim_state.trajectory_rotation = 0*[0 0 0; 20 20 20; 20 20 0; 10 20 0;  30 30 30; 0 0 0]*DEG2RAD;
sim_state.trajectory_accel = 0*[0 0 0; 0.25 0.25 0; 0.25 0.25 0; -0.25 -0.25 0; -2 -2 0; 0 0 0];
sim_state.velocity = zeros(sim_state.num_sim_samples, 3);
sim_state.position = zeros(sim_state.num_sim_samples, 3);
sim_state.accel = zeros(num_sim_samples, 3);
sim_state.orientation = zeros(num_sim_samples, 4);
sim_state.velocity(1,:) = initial_pos;
sim_state.position(1,:) = initial_vel;
sim_state.accel(1,:) = initial_accel;
sim_state.orientation(1,:) = eul2quat(sim_state.trajectory_rotation(1,1), sim_state.trajectory_rotation(1,2), sim_state.trajectory_rotation(1,3));

sim_state = init_simulator(sim_state, time_step, g);

f2 = figure; % for rotating cube
f3 = figure;

sensor_model = [];
sensor_model = init_sensor_model(sensor_model, time_step);

Qd =  [sensor_model.Q_q    zeros(3,3);
      zeros(3,3)          sensor_model.Q_gyro_bias];
      
P = [sensor_model.P_q    zeros(3,3);
	zeros(3,3)          sensor_model.P_gyro_bias];

ness = 0;

% initial estimates
sv = []; % state vector
sv.gyro_bias_est = 0*[0.0127 -0.0177 -0.0067]';
sv.accel_bias_est = [-0.06 0 0]';
sv.position_est = sim_state.position(1,:)'; % col vectors
sv.velocity_est = sim_state.velocity(1,:)';
q_prev = sim_state.orientation(1,:);
sv.q_est = q_prev';

accelUpdateFrequency = 10;

% visualization variables
position_history = [];
velocity_history = [];
roll_cov_history = [];
q_est_history = [];
roll_error = [];
pitch_error = [];
yaw_error = [];

% for recording video
% vidObj1 = VideoWriter('sf1_sim_rotating_cube.avi');
% open(vidObj1);
% vidObj2 = VideoWriter('sf1_sim_running_plot.avi');
% open(vidObj2);

for i = 2:num_sim_samples-500
	idx = i;
    % Get current sensor values from the simulator
	[omega_measured, accel_measured, q_prev, dcm_true, q_true] = sim_imu_tick(sim_state, time_step, idx, sensor_model, q_prev);
    
    % Calculate estimated values of omega and accel by subtracting the
    % estimated values of the biases
    omega_est = omega_measured - sv.gyro_bias_est;
    accel_est = accel_measured - sv.accel_bias_est;
 
    
    phi         = omega_est*time_step;
    delta_vel   = accel_est*time_step; % not necessary for accel-gyro sensor fusion
    phi3 = make_skew_symmetric_3(phi);
    
    % update current orientation estimate. We maintain orientation estimate
    % as a quaternion and obtain the corresponding DCM whenever needed
    sv.q_est = apply_small_rotation(phi, sv.q_est);
    sv.dcm_est = quat2dc(sv.q_est);
    
    % Integrate accel (after transforming to global frame) to obtain
    % velocity and position
    orig_velocity = sv.velocity_est;
    sv.velocity_est = sv.velocity_est + sv.dcm_est*delta_vel - g*time_step ;
    final_velocity = sv.velocity_est;
    sv.position_est = sv.position_est + ((orig_velocity + final_velocity)/2)*time_step;
    
    % State transition matrix
    F = eye(6) +   [-phi3       eye(3)*time_step;
                   zeros(3,6)];
 
    % propagate covariance 
    P = F*P*F' + Qd; 
  
    % For plotting purposes
    velocity_history(end+1,:) = sv.velocity_est';
    position_history(end+1,:) = sv.position_est';
    q_est_history(end+1, :) = sv.q_est';
    
    % Apply accelerometer measurements.
    if (mod(i, accelUpdateFrequency) == 0)
        % apply accel measurements. The updated state vector and covariance
        % matrix are returned. 
        [sv, P] = process_accel_measurement_update2(sv, accel_est, P, sensor_model.R_accel, g);
        gyro_bias_est 	= sv.gyro_bias_est;
		q_est           = sv.q_est;
        dcm_est         = sv.dcm_est;
    end
    
    % Rest is for visualization only
    roll_cov_history(end+1,:) = P(1,1);
    P_trace = sqrt(P(1,1)) %trace(P)/6;
    [r1 p1 y1] = dc2euler(dcm_true); % true euler angles
    [r2 p2 y2] = dc2euler(sv.dcm_est); % noisy, biased angles
    delta_rot = dcm_true*sv.dcm_est'-eye(3);
    roll_error(end+1) = norm([delta_rot(2,3)]); %yaw_error(end+1) = y2-y1; pitch_error(end+1) = p2-p1;
    
    % Visualization only: Shows a running plot of various states
    % show_running_plot(r1,p1,y1, r2, p2, y2, omega_measured, P_trace, f3);
    % Visualization only: Shows a totating cube corresponding to the
    % current estimated DCM
 %   draw_rotating_cube(sv.position_est, sv.dcm_est, f2); 
    % Visualization only:
    % write every third frame to make the video faster/smaller
%     if (mod(i, 3) == 0)        
%         writeVideo(vidObj1, getframe(f2));
%         writeVideo(vidObj2, getframe(f3));
%     end
    pause(0.01);
end

% close(vidObj1);
% close(vidObj2);
%figure;
%plot(position_history(:,1));
% solve for steady state P matrix
% H_accel = [make_skew_symmetric_3(g)*dcm_est zeros(3,3)];
% FF = F(7:12, 7:12);
% Qdd = Qd(7:12, 7:12);
% R = sensor_model.R_accel;
% [X,L,G] = dare(FF', H_accel', Qdd, R)
P
figure
plot(3*sqrt(roll_cov_history), 'r');
hold on;
plot(roll_error, 'b');
plot(-3*sqrt(roll_cov_history), 'r');
title({'Roll angle residual plot (Simulated)', 'Measurement Update Frequency = Sampling Frequency/10'}, 'interpreter', 'none', 'FontWeight','Normal');
ylabel('roll error sigma (radians)');
xlabel('time');
legend('3 sigma (roll)', 'roll angle residual', '-3 sigma (roll)');
ylim([-0.025, 0.025]);

var(roll_error)
mean(roll_error)
mean(yaw_error)
mean(pitch_error)
%save('simulated_imu_data', 'gyro_measurements', 'accel_measurements');
end


