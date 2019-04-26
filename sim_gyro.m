function sim_gyro()
close all;
time_step = 1/50;
time_to_run = 100; % seconds
num_points = round(time_to_run/time_step);
DEG2RAD = pi/180;
rotation_axis = [0 0 1];
% rotate at the rate of 5 degrees/sec
rot_delta = 5*time_step*DEG2RAD;
q_prev = [1 rotation_axis].*[1 0 0 0];
dir = 1;
q_est = q_prev';

global yaw_est_history;
global pitch_est_history;
global roll_est_history;

global yaw_true_history;
global pitch_true_history;
global roll_true_history;
global omega_measured_history;
global omega_true_history;
global state_uncertainity_history;
global NESS_history;

len = 100;

yaw_est_history = zeros(1,len);
pitch_est_history = zeros(1,len);
roll_est_history = zeros(1,len);
yaw_true_history = zeros(1,len);
pitch_true_history = zeros(1,len);
roll_true_history = zeros(1,len);
omega_measured_history = zeros(1,len);
omega_true_history = zeros(1,len);
state_uncertainity_history = zeros(1, len);
NESS_history = zeros(1, len);

f1 = figure; % for running plot
f2 = figure; % for rotating cube
f3 = figure;

trajectory = [0 0 0; 30 0 0; 30 30 0;  30 30 30; 0 0 0];
%trajectory = [0 0 0; 0 0 0];
trajectory = trajectory.*DEG2RAD;
% generate 20 sec worth of data at 100Hz.
%tq: trajectory quaternions
num_points = 2000;
num_wp = size(trajectory,1);
num_segments = num_wp - 1;
tq = [];
for i = 1: num_wp-1
    %begin quaternion
    qa = eul2quat(trajectory(i,1), trajectory(i,2), trajectory(i,3));
    qb = eul2quat(trajectory(i+1,1), trajectory(i+1,2), trajectory(i+1,3));
    num_points_segment = num_points/num_segments;
    for j = 1: num_points_segment
        tq(end+1,:) = slerp(qa, qb, j/num_points_segment);
    end
end

gyro_random_noise_sigma = [0.08 0.068 0.048]'; %rad/sec, representative values for MPU6050
gyro_bias_noise_sigma = 0.0001*ones(3,1); %rad/sec/sec
gyro_bias = [0.0127 0.0177 0.0067]'; % rad/sec 

% Process covariance
Q1 = diag([ 3.0983e-10   2.7251e-08     2.4144e-08]);
%Q1 = diag([0.08 0.068 0.048].*[0.08 0.068 0.048]);
Q = [Q1 zeros(3,3); zeros(3,3) 0.1*Q1];

% initial uncertainty
gyro_bias_est = [0 0 0]';
gyro_sigma = 0.1;
gyro_bias_sigma = 0.1;
P = [gyro_sigma^2*eye(3) zeros(3,3); zeros(3,3) gyro_bias_sigma^2*eye(3)];
P = zeros(6,6);
% measurement uncertainity
accel_noise_cov = [0.003 0.003 0.004]; % representative values for MPU6050
accel_noise_sigma = sqrt(accel_noise_cov)';
R = diag(accel_noise_cov);

q_prev = tq(1,:);
g = [0 0 1]';
err_state = zeros(6,1);
q_est = [1 0 0 0]';
roll_error = [];
gyro_measurements = [];
accel_measurements = [];
ness = 0;
phi_accum = [];
orientation_residual_accum = [];
accel_residual_accum = [];
for i = 2:num_points
    %rotation_axis = [0 0 1];
    % rotate at the rate of 5 degrees/sec
    %rot_delta = 5*time_step*DEG2RAD*dir;
    %sq = [1 rotation_axis].*[cos(rot_delta/2) sin(rot_delta/2) sin(rot_delta/2) sin(rot_delta/2)];
    %q = quatmult(q_prev, sq);
    %dq = (q-q_prev)./time_step;
    q = tq(i,:);
    dcm_true = quat2dc(q);
    accel = dcm_true'*g;
   % drawcube(dcm_true, f1);
    dq = (q - q_prev);
    q_prev = q;
    q_conj = [q(1) -q(2) -q(3) -q(4)];
    omega = 2*quatmult(q_conj, dq);
    omega = omega(2:4);
    % divide by dt to get angular rate
    omega = omega/time_step;
    % add noise
    gyro_bias_noisy = gyro_bias + normrnd(0, gyro_bias_noise_sigma);
    omega_measured = omega' + normrnd(0, gyro_random_noise_sigma);% + gyro_bias_noisy;
    accel_bias = [0 0 0];
    accel_measured = accel + normrnd(0, accel_noise_sigma) + accel_bias';
    gyro_measurements(end+1,:) = omega_measured';
    accel_measurements(end+1,:) = accel_measured';
    % Run KFplot
    omega_est = omega_measured - gyro_bias_est;
    %omega_est = omega_measured - gyro_bias;
    
    phi = omega_est*time_step;
    phi_accum(:,end+1) = omega_measured*time_step';
    phi3 = make_skew_symmetric_3(phi);
    phi4 = make_skew_symmetric_4(phi);
    s = norm(phi)/2;
    sin_sbys_approx = 1-s^2/factorial(3) + s^4/factorial(5) - s^6/factorial(7);
    exp_phi4 = eye(4)*cos(s) - 1/2*phi4*sin_sbys_approx;
    q_est = exp_phi4*q_est;
    F = eye(6) + [-phi3 eye(3)*time_step; zeros(3,3) zeros(3,3)];
    Qd = diag([0.08^2 0.068^2 0.0482^2 0.0001^2 0.0001^2 0.0001^2])*time_step;
    P1 = F*P*F' + Qd;    
    dcm_est = quat2dc([q_est(1), q_est(2), q_est(3), q_est(4)]);
    % calculate accel residual;
    % apply updates every 50 iterations to make it easy to see impact on
    % covariance
    % bias update
    if (i == num_points-1)
        k = 0;
    end 
    if (mod(i, 5) == 0)
        accel_residual = dcm_est*accel_measured - g;
        
        H_accel = [make_skew_symmetric_3(g)*dcm_est zeros(3,3)];
        K_gain = P1*H_accel'*inv([H_accel*P1*H_accel' + R]);
        xt = K_gain*accel_residual;
        P = P1 - K_gain*H_accel*P1;
        ness = xt(1:2)'*inv(P(1:2, 1:2))*xt(1:2)
        
        gyro_bias_est = gyro_bias_est - xt(4:6)
        % orientation update
        phi = [xt(1) xt(2) xt(3)];
        phi3 = make_skew_symmetric_3(phi);
        phi4 = make_skew_symmetric_4(phi);
        s = norm(phi)/2;
        sin_sbys_approx = 1-s^2/factorial(3) + s^4/factorial(5) - s^6/factorial(7);
        exp_phi4 = eye(4)*cos(s) - 1/2*phi4*sin_sbys_approx;
        q_est = exp_phi4*q_est;
    else
        P = P1;
    end
    P_trace = trace(P);
    %drawcube(dcm_est, f2);
    [r1 p1 y1] = quat2eul(q); % true euler angles
    [r2 p2 y2] = dc2euler(dcm_est); % noisy, biased angles
    show_running_plot(r1,p1,y1, r2, p2, y2, ness, omega_measured, omega, P_trace, len, f3);
    roll_error(end+1) = r2-r1;
    pause(0.01);
    %     omega_measured_integrated = omega_measured*time_step;
    %     phi = omega_measured_integrated;
    %     phi3 = make_skew_symmetric_3(phi);
    %     phi4 = make_skew_symmetric_4(phi);
    %     s = norm(phi)/2;
    %     sin_sbys_approx = 1-s^2/factorial(3) + s^4/factorial(5) - s^6/factorial(7);
    %     exp_phi4 = eye(4)*cos(s) - 1/2*phi4*sin_sbys_approx;
    %     q_est = exp_phi4*q_est;
    %     dcm_est = quat2dc(q_est);
    %     %drawcube(dcm_est, f2);
    %     [r1 p1 y1] = quat2eul(q); % true euler angles
    %     [r2 p2 y2] = dc2euler(dcm_est); % noisy, biased angles
    %     show_running_plot(r1,p1,y1, r2, p2, y2, omega_measured, omega, len, f1);
    %     pause(0.01);
end
save('simulated_imu_data', 'gyro_measurements', 'accel_measurements');
end

function q = quatmul(q1, q2)
q = [q1(4) q1(3) -q1(2) q1(1);
    -q1(3) q1(4) q1(1) q1(2);
    q1(2) -q1(1) q1(4) q1(3);
    -q1(1) -q1(2) -q1(3) q1(4)]*q2;
end

function drawcube(dcm, f)
% clear current plot
figure(f);
cla(gca);
% Define a six row by four column matrix to define the six cube faces
fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];

% Define an eight row by three column matrix to define the vertices at which
% the faces meet
vm = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];

% Plot the cube ----- gives each face a different color and creates the
% cube at a convenient viewing angle
clear cdata;
cdata = [
    0 0 0; % black
    1 0 0; % red
    1 0 1; % magenta
    0 0 1; % blue
    0 1 0; % green
    1 1 0; % yellow
    1 1 1; % white
    0 1 1; % cyan
    ];
vm = vm*dcm;
patch('Vertices',vm,'Faces',fm,'FaceVertexCData',cdata,'FaceColor','interp');

axis equal;
xlim([-2 2]); ylim([-2 2]); zlim([-2 2]);
view(3)
end

function show_running_plot(r1,p1,y1, r2, p2, y2, ness, omega_measured, omega_true, state_uncertainity, len, f)
global yaw_est_history;
global pitch_est_history;
global roll_est_history;
global yaw_true_history;
global pitch_true_history;
global roll_true_history;
global omega_measured_history;
global omega_true_history;
global state_uncertainity_history;
global NESS_history;

yaw_est_history(1:len-1) = yaw_est_history(2:len);
yaw_est_history(len) = y2;
pitch_est_history(1:len-1) = pitch_est_history(2:len);
pitch_est_history(len) = p2;

yaw_true_history(1:len-1) = yaw_true_history(2:len);
yaw_true_history(len) = y1;
pitch_true_history(1:len-1) = pitch_true_history(2:len);
pitch_true_history(len) = p1;

omega_measured_history(1:len-1) = omega_measured_history(2:len);
omega_measured_history(len) = omega_measured(3);

omega_true_history(1:len-1) = omega_true_history(2:len);
omega_true_history(len) = omega_true(3);

state_uncertainity_history(1:len-1) = state_uncertainity_history(2:len);
state_uncertainity_history(len) = state_uncertainity;

NESS_history(1:len-1) = NESS_history(2:len);
NESS_history(len) = ness;

figure(f);
subplot(4,1,1);
plot([1:len], yaw_est_history, 'b', [1:len], yaw_true_history, 'r');
axis([0 100 -1 1]);
subplot(4,1,2);
plot([1:len], pitch_est_history , 'b', [1:len], pitch_true_history, 'r');
axis([0 100 -1 1]);
subplot(4,1,3);
plot([1:len], omega_measured_history , 'b', [1:len], omega_true_history, 'r');
axis([0 100 -0.1 0.1]);
subplot(4,1,4);
plot([1:len], state_uncertainity_history , 'b');
axis([0 100 0 0.5]);
end

function qm = slerp(qa, qb, t)
% quaternion to return
% Calculate angle between them.
cosHalfTheta = qa(1) * qb(1) + qa(2) * qb(2) + qa(3) * qb(3) + qa(4) * qb(4);
%if qa=qb or qa=-qb then theta = 0 and we can return qa
if (abs(cosHalfTheta) >= 1.0)
    qm(1) = qa(1);
    qm(2) = qa(2); qm(3) = qa(3); qm(4) = qa(4);
    return;
end
% Calculate temporary values.
halfTheta = acos(cosHalfTheta);
sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta);
% if theta = 180 degrees then result is not fully defined
%we could rotate around any axis normal to qa or qb
if (abs(sinHalfTheta) < 0.001)% fabs is floating point absolute
    qm(1) = (qa(1) * 0.5 + qb(1) * 0.5);
    qm(2) = (qa(2) * 0.5 + qb(2) * 0.5);
    qm(3) = (qa(3) * 0.5 + qb(3) * 0.5);
    qm(4) = (qa(4) * 0.5 + qb(4) * 0.5);
    return;
end
ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
ratioB = sin(t * halfTheta) / sinHalfTheta;
% calculate Quaternion.
qm(1) = (qa(1) * ratioA + qb(1) * ratioB);
qm(2) = (qa(2) * ratioA + qb(2) * ratioB);
qm(3) = (qa(3) * ratioA + qb(3) * ratioB);
qm(4) = (qa(4) * ratioA + qb(4) * ratioB);
end

function phi4 = make_skew_symmetric_4(phi)
phi4 = [0       phi(1)      phi(2)      phi(3);
    -phi(1)     0           -phi(3)     phi(2);
    -phi(2)     phi(3)      0           -phi(1);
    -phi(3)     -phi(2)     phi(1)      0];
end

function phi3 = make_skew_symmetric_3(phi)
phi3 = [0       -phi(3)      phi(2);
    phi(3)     0           -phi(1);
    -phi(2)     phi(1)      0];
end