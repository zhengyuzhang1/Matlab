sz = [1 0;0 -1];
sx = [0 1;1 0];
sy = [0 -1i;1i 0];

w00 = 505*2*pi*0.00107e6;
w1 = 2 * pi * 0.13e6;
A = 2 * pi * 13.7e6;
th = 7/180*pi;

Az = A * cos(th);
Ax = A * sin(th);
w01 = sqrt((w00 + Az)^2 + Ax^2);
sth = Ax/w01;
cth = (w00+Az)/w01;
szp = cth*sz+sth*sx;

H0 = w00 * sz /2;
H1 = w01 * (sth * sx + cth * sz) / 2;

n1 = [1 0 0];
t = 1e-7:1e-7:1e-6;
trfovert = 0.1:0.1:1;%targetangle / Rabi_f / 2;
w = linspace(w00,w01,10);

a0 = zeros(length(t),length(trfovert),length(w),3);
a1 = zeros(size(a0));
ax = zeros(length(t),length(trfovert),length(w));
angle0 = zeros(size(ax));
angle1 = zeros(size(ax));
phi0 = 0;
n00 = [0 0 1];
n01 = [sth 0 cth];
for i = 1:length(t)
    for j = 1:length(trfovert)
        trf = trfovert(j)*t(i);
        for k = 1:length(w)
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
        end
    end
end
%%
plotx = zeros(length(t)*length(trfovert)*length(w),1);
ploty = zeros(size(plotx));
plotz = zeros(size(plotx));
plotc = zeros(size(plotx));
ind = 0;
for i = 1:length(t)
    for j = 1:length(trfovert)
        for k = 1:length(w)
            ind = ind + 1;
            plotx(ind) = t(i);
            ploty(ind) = trfovert(j)*t(i);
            plotz(ind) = w(k);
            plotc(ind) = ax(i,j,k);
        end
    end
end
%figure;scatter(ploty,plotz,5,plotc);
figure;scatter3(plotx,ploty,plotz,5,plotc);