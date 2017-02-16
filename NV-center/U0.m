function f = U0(w, azx, t)
f = expm(-0.5i*w*t*[1,0;0,-1])*expm(-1i*t*((azx(1)+w)*[1,0;0,-1]+azx(2)*[0,1;1,0]))*expm(-0.5i*w*t*[1,0;0,-1]);
end