function U1 = U1rf(Az, Ax, w, rabi, w_rf, t_rf, t, b)

bx = b(1) / norm(b);
by = b(2) / norm(b);
bz = b(3) / norm(b);
sz = [1 0;0 -1];
sx = [0 1;1 0];
sy = [0 -1i;1i 0];
H1 = (w+Az)/2*sz+Ax/2*sx;
w1 = sqrt((w+Az)^2+Ax^2);
sth = Ax/w1;
cth = (w+Az)/w1;
szp = cth*sz+sth*sx;
sxp = (cth^2*bx-cth*sth*bz)*sx+by*sy+(sth^2*bz-cth*sth*bx)*sz;
U1 = expm(-1i*(t-t_rf)*H1)*expm(-1i*t_rf*w_rf/2*szp)...
    *expm(-1i*t_rf*((w1-w_rf)/2*szp+rabi/2*sxp));

end    