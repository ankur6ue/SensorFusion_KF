function qm = interpolate_quat(qa, qb, t)
% interpolating from quaternion qa to qb at time t using slerp    
qm = slerp(qa, qb, t);
end
function qm = slerp(qa, qb, t)
% quaternion to return
% Calculate angle between them.
cosHalfTheta = qa(1) * qb(1) + qa(2) * qb(2) + qa(3) * qb(3) + qa(4) * qb(4);
%if qa=qb or qa=-qb then theta = 0 and we can return qa
if (abs(cosHalfTheta) >= 1.0)
    qm(1) = qa(1);
    qm(2) = qa(2); qm(3) = qa(3); qm(4) = qa(4);
    return;
end
% Calculate temporary values.
halfTheta = acos(cosHalfTheta);
sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta);
% if theta = 180 degrees then result is not fully defined
%we could rotate around any axis normal to qa or qb
if (abs(sinHalfTheta) < 0.001)% fabs is floating point absolute
    qm(1) = (qa(1) * 0.5 + qb(1) * 0.5);
    qm(2) = (qa(2) * 0.5 + qb(2) * 0.5);
    qm(3) = (qa(3) * 0.5 + qb(3) * 0.5);
    qm(4) = (qa(4) * 0.5 + qb(4) * 0.5);
    return;
end
ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
ratioB = sin(t * halfTheta) / sinHalfTheta;
% calculate Quaternion.
qm(1) = (qa(1) * ratioA + qb(1) * ratioB);
qm(2) = (qa(2) * ratioA + qb(2) * ratioB);
qm(3) = (qa(3) * ratioA + qb(3) * ratioB);
qm(4) = (qa(4) * ratioA + qb(4) * ratioB);
end
