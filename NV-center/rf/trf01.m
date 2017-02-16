function U = trf01(w_rf,Az,Ax,w,t)
% transforming from rotating frame, w_rf*sz to w_rf*(cos(th)*sz+sin(th)*sx),each
% rotated total time t

sz = [1 0;0 -1];
sx = [0 1;1 0];
th = atan(Ax / (Az + w));
U = expm(1i*t*w_rf/2*(cos(th)*sz+sin(th)*sx))*expm(-1i*t*w_rf/2*sz);

end