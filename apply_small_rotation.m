function q_est = apply_small_rotation(phi, q_est)
s = norm(phi)/2;
x = 2*s;
phi4 = make_skew_symmetric_4(phi);
% approximation for sin(s)/s
sin_sbys_approx = 1-s^2/factorial(3) + s^4/factorial(5) - s^6/factorial(7);
exp_phi4 = eye(4)*cos(s) + 1/2*phi4*sin_sbys_approx;
% method 1
q_est1 = exp_phi4*q_est;
% verify: should be equal to:
% method 2 (direct matrix exponential)
q_est2 = expm(0.5*(phi4))*q_est;
% verify should be equal to:
% method 3: propagate by quaternion multiplication (doesn't preserve unit
% norm)
dq = [1 phi(1)/2 phi(2)/2 phi(3)/2];
dq_conj = [1 -phi(1)/2 -phi(2)/2 -phi(3)/2]';
q_est3 = quatmul(q_est, dq'); 
q_est = q_est1;
end