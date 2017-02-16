function dy = rfH(t, y, w, rabi, wrf, phi)
H = w / 2 * [1,0;0,-1] + (rabi' * cos(wrf * t + phi)) * [0,1;1,0];
dy = -1i * H * y;
end