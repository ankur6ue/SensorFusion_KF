function sv = update_state(sv, time_step, g, phi, vel)
sv.q_est = apply_small_rotation(phi, sv.q_est);
sv.dcm_est = quat2dc(sv.q_est);
orig_velocity = sv.velocity_est;
sv.velocity_est = sv.velocity_est + sv.dcm_est*vel - g*time_step ;
final_velocity = sv.velocity_est;
sv.position_est = sv.position_est + ((orig_velocity + final_velocity)/2)*time_step;
sv.p_C_in_G_est = sv.position_est;
end