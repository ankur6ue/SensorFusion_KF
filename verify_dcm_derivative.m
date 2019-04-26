% Original vector v
v = [1 1 1];
% Convert to unit vector
v = v./norm(v);
% Initial rotation (euler angles)
eul = [30 30 30]/57.3;
% Get quaternion and DCM representation for this rotation
q0 = eul2quat(30/57.3, 30/57.3, 30/57.3);
C0 = euler2dc(30/57.3, 30/57.3, 30/57.3);
% Apply rotation
v1 = C0*v';
% Small rotation (in radians)
phi = [0.01 0 0.01];
% Obtain the 3*3 and 4*4 skew symmetric matrices
phi3 = make_skew_symmetric_3(phi);
phi4 = make_skew_symmetric_4(phi);
s = norm(phi)/2;
x = norm(phi);
exp_phi4 = eye(4)*cos(s) - 1/2*phi4*sin(s)/s;
% Propagage quaternion
q1 = exp_phi4*q0';
% Propagate DCM (Rodrigues Formula)
C1 = C0*(eye(3) + sin(x)/x*phi3 + (1-cos(x))/x^2*phi3*phi3);
% Obtain the DCM corresponding to q1 to compare with C1
C_q1 = quat2dc(q1);
% Verify C_q1 ~ C1
% Now verify derivative wrt phi
% True value of the derivative
d1 = C1*v'-C0*v'
d2 = C0*make_skew_symmetric_3(v')*phi'
% verify d1 ~ d2

% d1 should equal -
function phi4 = make_skew_symmetric_4(phi)
phi4 = [0       phi(1)      phi(2)      phi(3);
    -phi(1)     0           -phi(3)     phi(2);
    -phi(2)     phi(3)      0           -phi(1);
    -phi(3)     -phi(2)     phi(1)      0];
end

function phi3 = make_skew_symmetric_3(phi)
phi3 = [0       -phi(3)      phi(2);
    phi(3)     0           -phi(1);
    -phi(2)     phi(1)      0];
end