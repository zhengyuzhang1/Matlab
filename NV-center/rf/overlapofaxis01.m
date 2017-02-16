function f = overlapofaxis01(t, trf, w00, w01, w1, w, phi, n0, n1, extraU)

v01 = expm(-1i*(t(i)-trf)*H0)*anyrf(trf,w00,w1,w(k),phi0,n00,n1);
v02 = expm(-1i*(t(i)-trf)*H1)*anyrf(trf,w01,w1,w(k),phi0+w(k)*t(i),n01,n1);
v03 = expm(-1i*(t(i)-trf)*H0)*anyrf(trf,w00,w1,w(k),phi0+w(k)*2*t(i),n00,n1);
v04 = expm(-1i*(t(i)-trf)*H1)*anyrf(trf,w01,w1,w(k),phi0+w(k)*3*t(i),n01,n1);
u0 = v04 * v03 * v02 * v01;
%u0 = expm(1i * 4*t(i) * w01*(sth*sx+cth*sz)/2)*u0;
v11 = expm(-1i*(t(i)-trf)*H1)*anyrf(trf,w01,w1,w(k),phi0,n01,n1);
v12 = expm(-1i*(t(i)-trf)*H0)*anyrf(trf,w00,w1,w(k),phi0+w(k)*t(i),n00,n1);
v13 = expm(-1i*(t(i)-trf)*H1)*anyrf(trf,w01,w1,w(k),phi0+w(k)*2*t(i),n01,n1);
v14 = expm(-1i*(t(i)-trf)*H0)*anyrf(trf,w00,w1,w(k),phi0+w(k)*3*t(i),n00,n1);
u1 = v14 * v13 * v12 * v11;
%u1 = expm(1i * 4*t(i) * w01*(sth*sx+cth*sz)/2)*u1;
[a0(i,j,k,:),angle0(i,j,k)] = axisangle(u0);
[a1(i,j,k,:),angle1(i,j,k)] = axisangle(u1);
ax(i,j,k) = reshape(a0(i,j,k,:),1,3)* reshape(a1(i,j,k,:),3,1);