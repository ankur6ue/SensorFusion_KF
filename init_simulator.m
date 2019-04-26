function [sim_state] = init_simulator(sim_state, time_step, g)
num_wp = size(sim_state.trajectory_rotation,1);
num_segments = num_wp - 1;
tq = [];
a_prev = 0;
ap = 0;
idx = 2;
for i = 1: num_wp-1
    % beginning quaternion
    qa = eul2quat(sim_state.trajectory_rotation(i,1), sim_state.trajectory_rotation(i,2), sim_state.trajectory_rotation(i,3));
    % ending quaternion
    qb = eul2quat(sim_state.trajectory_rotation(i+1,1), sim_state.trajectory_rotation(i+1,2), sim_state.trajectory_rotation(i+1,3));
    ab = sim_state.trajectory_accel(i,:); % sim_state.accel beginning
    ae = sim_state.trajectory_accel(i+1,:); % sim_state.accel end
    % number of points on each segment of the trajectory
    num_points_segment = sim_state.num_sim_samples/num_segments;
    for j = 1: num_points_segment
        sim_state.orientation(idx,:) = interpolate_quat(qa, qb, j/num_points_segment);
        ap = ap + (ae-ab)/num_points_segment;
        sim_state.accel(idx,:) = ap + g';
        sim_state.velocity(idx,:) = sim_state.velocity(idx-1,:)+(a_prev+ap)/2*time_step;
        sim_state.position(idx,:) = sim_state.position(idx-1,:) + (sim_state.velocity(idx,:)+sim_state.velocity(idx-1,:))/2*time_step;
        a_prev = ap;
		idx = idx + 1;
    end
end
end
