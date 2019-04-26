
function show_running_plot(r1,p1,y1, r2, p2, y2, omega_measured, state_uncertainity, f)
global state_history;
len = state_history.len;

state_history.yaw_est(1:len-1) = state_history.yaw_est(2:len);
state_history.yaw_est(len) = y2;
state_history.pitch_est(1:len-1) = state_history.pitch_est(2:len);
state_history.pitch_est(len) = p2;
state_history.roll_est(1:len-1) = state_history.roll_est(2:len);
state_history.roll_est(len) = r2;

state_history.yaw_true(1:len-1) = state_history.yaw_true(2:len);
state_history.yaw_true(len) = y1;
state_history.pitch_true(1:len-1) = state_history.pitch_true(2:len);
state_history.pitch_true(len) = p1;
state_history.roll_true(1:len-1) = state_history.roll_true(2:len);
state_history.roll_true(len) = r1;

state_history.omega_measured(1:len-1) = state_history.omega_measured(2:len);
state_history.omega_measured(len) = omega_measured(3);

state_history.state_uncertainity(1:len-1) = state_history.state_uncertainity(2:len);
state_history.state_uncertainity(len) = state_uncertainity;

figure(f);
subplot(4,1,1);
plot([1:len], state_history.yaw_est - state_history.yaw_true, 'b');
title('yaw error (radians)','FontWeight','Normal');
axis([0 100 -0.2 0.2]);
subplot(4,1,2);
plot([1:len], state_history.pitch_est - state_history.pitch_true , 'b');
title('pitch error (radians)', 'FontWeight','Normal');
axis([0 100 -0.2 0.2]);
subplot(4,1,3);
plot([1:len], state_history.roll_est - state_history.roll_true , 'b');
title('roll error (radians)', 'FontWeight','Normal');
axis([0 100 -0.2 0.2]);
subplot(4,1,4);
plot([1:len], state_history.state_uncertainity , 'b');
title('State Covariance (trace(P))', 'FontWeight','Normal');
axis([0 100 0 0.1]);

end
