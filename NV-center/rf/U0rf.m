function U0 = U0rf(w, rabi, w_rf, t_rf, t, b)
%in lab frame

bx = b(1) / norm(b);
by = b(2) / norm(b);
sz = [1 0;0 -1];
sx = [0 1;1 0];
sy = [0 -1i;1i 0];
U0 = expm(-1i*((t-t_rf)*w+t_rf*w_rf)/2*sz)*expm(-1i*t_rf*((w-w_rf)/2*sz+rabi/2*(bx*sx+by*sy)));

end
    