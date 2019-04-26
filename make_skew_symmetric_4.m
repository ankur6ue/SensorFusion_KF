
function phi4 = make_skew_symmetric_4(phi)
phi4 = [0       -phi(1)      -phi(2)      -phi(3);
    phi(1)     0           phi(3)     -phi(2);
    phi(2)     -phi(3)      0           phi(1);
    phi(3)     phi(2)     -phi(1)      0];
end
