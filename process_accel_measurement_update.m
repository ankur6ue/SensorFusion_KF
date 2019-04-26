function [sv, P] = process_accel_measurement_update(sv, accel_est, P, R_accel, g)
q_est   		= sv.q_est;
dcm_est 		= sv.dcm_est;
gyro_bias_est 	= sv.gyro_bias_est;

accel_residual = dcm_est*accel_est - g;
H_accel = [zeros(3,3) zeros(3,3) make_skew_symmetric_3(g)*dcm_est zeros(3,3) zeros(3,3)];
K_gain = P*H_accel'*inv([H_accel*P*H_accel' + R_accel]);
x_corr = K_gain*accel_residual;
P = P - K_gain*H_accel*P;
%ness = x_corr(1:2)'*inv(P(1:2, 1:2))*x_corr(1:2);
x_est = apply_accel_update(x_corr(7:12), [q_est; gyro_bias_est]);
gyro_bias_est(1:2) = x_est(5:6);
q_est = x_est(1:4);

sv.q_est         = q_est;
sv.dcm_est       = dcm_est;
sv.gyro_bias_est = gyro_bias_est;

end


function x_est = apply_accel_update(x_corr, x_est)
% gyro bias
x_est(5:7) = x_est(5:7) - x_corr(4:6);
% orientation update
phi = [x_corr(1) x_corr(2) x_corr(3)];
phi3 = make_skew_symmetric_3(phi);
phi4 = make_skew_symmetric_4(phi);
s = norm(phi)/2;
sin_sbys_approx = 1-s^2/factorial(3) + s^4/factorial(5) - s^6/factorial(7);
exp_phi4 = eye(4)*cos(s) - 1/2*phi4*sin_sbys_approx;
x_est(1:4) = exp_phi4*x_est(1:4);
end
