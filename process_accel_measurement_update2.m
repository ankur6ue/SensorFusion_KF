% reduced model, applies to accel+gyro only. No camera
function [state, P] = process_accel_measurement_update2(state, accel_est, P, R_accel, g)
% current estimated states
q_est   = state.q_est;
dcm_est = state.dcm_est;
gyro_bias_est = state.gyro_bias_est;
% Calculate residual (in the global reference frame)
accel_residual = dcm_est*accel_est - g;
% Construct the H matrix
H_accel = [-make_skew_symmetric_3(g)*dcm_est zeros(3,3)];
% Calculate Kalman Gain
K_gain = P*H_accel'*inv([H_accel*P*H_accel' + R_accel]);
% Calculate the correction to the error state
x_corr = K_gain*accel_residual;
% Calculate the new covariance matrix
P = P - K_gain*H_accel*P;
% Apply the updates to get the new state
x_est = apply_accel_update(x_corr, [q_est; gyro_bias_est]);

gyro_bias_est(1:3) = x_est(5:7);
q_est = x_est(1:4);
%ness = [x_corr(1:2) x_corr(4:5)]'*inv(P(1:2, 1:2))*x_corr(1:2);
ness = x_corr'*inv(P)*x_corr;

% Update state vector
state.q_est         = q_est;
state.dcm_est       = dcm_est;
state.gyro_bias_est = gyro_bias_est;
end

function x_est = apply_accel_update(x_corr, x_est)
% gyro bias
x_est(5:7) = x_est(5:7) + x_corr(4:6);
% orientation update
phi =  -[x_corr(1) x_corr(2) x_corr(3)];
phi3 = make_skew_symmetric_3(phi);
phi4 = make_skew_symmetric_4(phi);
s = norm(phi)/2;
sin_sbys_approx = 1-s^2/factorial(3) + s^4/factorial(5) - s^6/factorial(7);
exp_phi4 = eye(4)*cos(s) + 1/2*phi4*sin_sbys_approx;
x_est(1:4) = exp_phi4*x_est(1:4);
end
