function q = quatmul(q1, q2)
% for the case when quaternion representation is [vector scalar]
% q = [q1(4) q1(3) -q1(2) q1(1);
%     -q1(3) q1(4) q1(1) q1(2);
%     q1(2) -q1(1) q1(4) q1(3);
%     -q1(1) -q1(2) -q1(3) q1(4)]*q2;
% for the case when quaternion representation is [scalar vector]
q = [q1(1) -q1(2) -q1(3) -q1(4);
    q1(2) q1(1) -q1(4) q1(3);
    q1(3) q1(4) q1(1) -q1(2);
    q1(4) -q1(3) q1(2) q1(1)]*q2;

end