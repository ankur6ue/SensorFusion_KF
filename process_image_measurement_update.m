function  [sv, P] =  process_image_measurement_update(sv, p_f_in_G, measurements, numFeatures, P, K, imageW, imageH, sigma)
p_C_in_G_est 	= sv.p_C_in_G_est;
velocity_est 	= sv.velocity_est;
gyro_bias_est 	= sv.gyro_bias_est;
accel_bias_est 	= sv.accel_bias_est;
q_est 			= sv.q_est;
dcm_est 		= sv.dcm_est;
position_est 	= sv.position_est;
KF_SV_Offset 	= sv.KF_SV_Offset;

global f4;
numVisible = 0;
H_image = [];
counter = 0;
residual = [];
visible = [];
for fid = 1:numFeatures
    if (measurements(2*fid-1) > 0 && measurements(2*fid-1) < imageW ...
            && measurements(2*fid) > 0 && measurements(2*fid) < imageH )
        visible(end+1) = fid;
        counter = counter + 1;
        numVisible = numVisible + 1;        
    end
end
% Only apply image measurement update if at least two feature points are visible. 
if (numVisible >= 2)
    % Measurement covariance matrix
    R_im = sigma^2*eye(numVisible*2);
    R_inv = inv(R_im);
    CF = 500;
    numIter = 0;
    H_image = [];
    residual = [];
    % Iterate until the CF < 0.25 or max number of iterations is reached
    while(CF > 0.05 && numIter < 10)
        numIter = numIter + 1;
        estimated_visible = [];
        measurements_visible = [];
        residual = [];
        for indx = 1:numVisible
            % index of this feature point
            fid = visible(indx);
            % estimated position in the camera coordinate system
            fi_e = dcm_est'*(p_f_in_G(fid,:)' - p_C_in_G_est); % feature in image, estimated
			% focal lengths (recall that the x axis of our coordinate
			% system points along the optical axis)
            fy = K(1,1); fz = K(2,2);
            J = 	[-fy*fi_e(2)/fi_e(1).^2  fy*1/fi_e(1)       0;
                    -fz*fi_e(3)/fi_e(1).^2    	0               fz*1/fi_e(1)];
            
            % Measurement matrix
            H_image(2*indx-1: 2*indx,:) = [-J*dcm_est' zeros(2,3) ...
                J*dcm_est'*make_skew_symmetric_3((p_f_in_G(fid,:)' - p_C_in_G_est)) zeros(2,3) zeros(2,3)];
            % estimated image measurement
            estimated_visible(2*indx-1: 2*indx) = K*[fi_e(2)/fi_e(1) fi_e(3)/fi_e(1) 1]';
            % actual measurement
            measurements_visible(2*indx-1: 2*indx) = measurements(2*fid-1: 2*fid)';
        end
        % vector of residuals
        residual = estimated_visible - measurements_visible;
        %show_image_residuals(f4, measurements_visible, estimated_visible);
        % Kalman gain
        K_gain = P*H_image'*inv([H_image*P*H_image' + R_im]);
        % Correction vector
        x_corr = K_gain*residual';
        % Updated covariance
        P = P - K_gain*H_image*P;
        % Apply image update and correct the current estimates of position,
        % velocity, orientation, gyro/accel bias
        x_est = apply_image_update(KF_SV_Offset, x_corr, ...
            [position_est; velocity_est; q_est; gyro_bias_est; accel_bias_est]);
        
        position_est    = x_est(KF_SV_Offset.pos_index:KF_SV_Offset.pos_index+KF_SV_Offset.pos_length-1);
        velocity_est    = x_est(KF_SV_Offset.vel_index:KF_SV_Offset.vel_index+KF_SV_Offset.vel_length-1);
        gyro_bias_est   = x_est(KF_SV_Offset.gyro_bias_index:KF_SV_Offset.gyro_bias_index + KF_SV_Offset.gyro_bias_length-1);
        accel_bias_est  = x_est(KF_SV_Offset.accel_bias_index:KF_SV_Offset.accel_bias_index + KF_SV_Offset.accel_bias_length-1);
        q_est           = x_est(KF_SV_Offset.orientation_index:KF_SV_Offset.orientation_index+KF_SV_Offset.orientation_length-1);

        dcm_est = quat2dc(q_est);
        p_C_in_G_est = position_est;
        
        % Cost function used to end the iteration
        CF = x_corr'*inv(P)*x_corr + residual*R_inv*residual';
    end
end

sv.p_C_in_G_est = p_C_in_G_est; 
sv.velocity_est = velocity_est;
sv.position_est = position_est;
sv.gyro_bias_est = gyro_bias_est;
sv.accel_bias_est = accel_bias_est;
sv.q_est = q_est;
sv.dcm_est = dcm_est;
end
