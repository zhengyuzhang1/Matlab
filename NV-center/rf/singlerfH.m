function dy = singlerfH(t, y, w, rabi)
H = w / 2 * [1,0;0,-1] + rabi * cos(w * t) * [0,1;1,0];
dy = -1i * H * y;
end