function image_measurements = sim_camera_tick(dcm_true, p_f_in_G, p_C_in_G, K, numFeatures, sigma_image)
for fid = 1: numFeatures
    % feature position in camera frame
    fic = dcm_true'*(p_f_in_G(fid,:) - p_C_in_G')'; 
    % apply perspective projection and add noise
    image_measurements(2*fid-1:2*fid) = K*[fic(2)/fic(1) fic(3)/fic(1) 1]' ...
        + normrnd(0, [sigma_image sigma_image]');
end
end