function sensor_model = init_sensor_model(sensor_model, time_step)
sensor_model.accel_noise_sigma = [0.003 0.003 0.004]'; % representative values for MPU6050
sensor_model.accel_noise_cov = (sensor_model.accel_noise_sigma).^2;
sensor_model.accel_bias = [-0.06 0 0]'; % actual MPU6050 bias

sensor_model.accel_bias_noise_sigma = [0.001 0.001 0.001];
sensor_model.R_accel = diag(sensor_model.accel_noise_cov);
sensor_model.Q_accel_bias = diag(sensor_model.accel_bias_noise_sigma.^2);
sensor_model.P_accel_bias = 0.1*eye(3,3);

sensor_model.gyro_random_noise_sigma = [5.4732e-04 6.1791e-04 6.2090e-04]'; %rad/sec, representative values for MPU6050
sensor_model.gyro_bias_noise_sigma = 0.00001*ones(3,1); %rad/sec/sec. Any small value appears to work
sensor_model.gyro_bias = [0.0127 -0.0177 -0.0067]'; % rad/sec
sensor_model.Q_gyro_bias = diag(sensor_model.gyro_bias_noise_sigma.^2);

% position uncertainity 
sensor_model.Q_p = diag(sensor_model.accel_noise_cov)*time_step;
% velocity uncertainity
sensor_model.Q_v = diag(sensor_model.accel_noise_cov)*time_step;
% quaternion (orientation) uncertainity
sensor_model.Q_q = diag(sensor_model.gyro_random_noise_sigma.^2)*time_step;

% Initial values for various elements of the covariance matrix. Only
% effects the transient behaviour of the filter, not the steady state
% behaviour. 

sensor_model.P_p = 0.1*eye(3,3);

sensor_model.P_v = 0.1*eye(3,3);

sensor_model.P_q = 0.1*eye(3,3);

sensor_model.P_gyro_bias = 0.1*eye(3,3);

sensor_model.P_accel_bias = 0.1*eye(3,3);

end