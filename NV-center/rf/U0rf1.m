function U0 = U0rf1(w, rabi, w_rf, t)

sz = [1 0;0 -1];
sx = [0 1;1 0];

U0 = expm(-1i*t*((w-w_rf)/2*sz+ rabi/2*sx));

end
    