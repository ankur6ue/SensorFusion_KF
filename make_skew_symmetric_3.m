
function phi3 = make_skew_symmetric_3(phi)
phi3 = [0       -phi(3)      phi(2);
    phi(3)     0           -phi(1);
    -phi(2)     phi(1)      0];
end