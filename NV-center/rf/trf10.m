function U = trf10(w_rf,Az,Ax,w,t)
% transforming from rotating frame, w_rf*(cos(th)*sz+sin(th)*sx) to w_rf*sz,each
% rotated total time t

sz = [1 0;0 -1];
sx = [0 1;1 0];
th = atan(Ax / (Az + w));
U = expm(1i*t*w_rf/2*sz)*expm(-1i*t*w_rf/2*(cos(th)*sz+sin(th)*sx));

end