
function x_est = apply_image_update(KF_SV_Offset, x_corr, x_est)

% update position
x_est(KF_SV_Offset.pos_index: KF_SV_Offset.pos_index+KF_SV_Offset.pos_length-1) = x_est(KF_SV_Offset.pos_index: KF_SV_Offset.pos_index+KF_SV_Offset.pos_length-1)...
    - x_corr(KF_SV_Offset.pos_index: KF_SV_Offset.pos_index+KF_SV_Offset.pos_length-1);

% update veloctiy
x_est(KF_SV_Offset.vel_index: KF_SV_Offset.vel_index+KF_SV_Offset.vel_length-1) = x_est(KF_SV_Offset.vel_index: KF_SV_Offset.vel_index+KF_SV_Offset.vel_length-1)...
    - x_corr(KF_SV_Offset.vel_index: KF_SV_Offset.vel_index+KF_SV_Offset.vel_length-1);

% orientation update
phi = -[x_corr(KF_SV_Offset.orientation_index) x_corr(KF_SV_Offset.orientation_index+1) x_corr(KF_SV_Offset.orientation_index+2)];

x_est(KF_SV_Offset.orientation_index:KF_SV_Offset.orientation_index+KF_SV_Offset.orientation_length-1) = ...
apply_small_rotation(phi, x_est(KF_SV_Offset.orientation_index:KF_SV_Offset.orientation_index+KF_SV_Offset.orientation_length-1));

% gyro bias update. -1 has to be added to the gyro/accel bias index because the length of the orientation
% correction is 3 while the length of the corresponding orientation estimate is 4 (quaternion)
x_est(KF_SV_Offset.gyro_bias_index: KF_SV_Offset.gyro_bias_index+KF_SV_Offset.gyro_bias_length-1) = x_est(KF_SV_Offset.gyro_bias_index: KF_SV_Offset.gyro_bias_index+KF_SV_Offset.gyro_bias_length-1)...
    + x_corr(KF_SV_Offset.gyro_bias_index-1: KF_SV_Offset.gyro_bias_index-1+KF_SV_Offset.gyro_bias_length-1);

% accel bias update
x_est(KF_SV_Offset.accel_bias_index: KF_SV_Offset.accel_bias_index+KF_SV_Offset.accel_bias_length-1) = x_est(KF_SV_Offset.accel_bias_index: KF_SV_Offset.accel_bias_index+KF_SV_Offset.accel_bias_length-1)...
    - x_corr(KF_SV_Offset.accel_bias_index-1: KF_SV_Offset.accel_bias_index-1+KF_SV_Offset.accel_bias_length-1);

end
