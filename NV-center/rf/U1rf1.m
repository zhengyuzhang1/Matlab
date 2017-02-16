function U1 = U1rf1(Az, Ax, w, rabi, w_rf, t)

sz = [1 0;0 -1];
sx = [0 1;1 0];
w1 = sqrt((w+Az)^2+Ax^2);
sintheta = Ax/w1;
costheta = (w+Az)/w1;
szp = costheta*sz+sintheta*sx;
sxp = -sintheta*sz+costheta*sx;
rabi_v = rabi*costheta;
U1 = expm(-1i*t*((w1-w_rf)/2*szp+rabi_v/2*sxp));

end
    