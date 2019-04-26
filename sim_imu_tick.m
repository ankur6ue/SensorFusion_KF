function [omega_measured, accel_measured, q_prev, dcm_true, q] = sim_imu_tick(sim_state, time_step, idx, sensor_model, q_prev)
if (idx < sim_state.num_sim_samples)
    q = sim_state.orientation(idx,:);
    dcm_true = quat2dc(q);
    % dcm rotates from IMU to Global frame, dcm' does the opposite
    dq = (q - q_prev);
    q_prev = q;
    q_conj = [q(1) -q(2) -q(3) -q(4)];
    % calculate the angular velocity required to go from one orientation to the next 
    omega = 2*quatmult(q_conj, dq)/time_step;
    % take the vector part
    omega = omega(2:4);
    % add noise
    gyro_bias_noisy = sensor_model.gyro_bias + normrnd(0, sensor_model.gyro_bias_noise_sigma);
    omega_measured = omega' + normrnd(0, sensor_model.gyro_random_noise_sigma)+ gyro_bias_noisy;
    accel_measured = dcm_true'*sim_state.accel(idx,:)' + sensor_model.accel_bias + normrnd(0, sensor_model.accel_noise_sigma);
    %   gyro_measurements(end+1,:) = omega_measured';
    %   accel_measurements(end+1,:) = accel_measured';
end
end