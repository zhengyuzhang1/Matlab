sz = [1 0;0 -1];
sx = [0 1;1 0];
sy = [0 -1i;1i 0];

w = 505*2*pi*0.00107e6;
rabi = 2 * pi * 0.13e6;
A = 2 * pi * 13.7e6;
Az = A * cos(7/180*pi);
Ax = A * sin(7/180*pi);
w1 = sqrt((w + Az)^2 + Ax^2);
sth = Ax/w1;
cth = (w+Az)/w1;
szp = cth*sz+sth*sx;
% H0 = 0.5 * w * sz;
% H1 = 0.5 * (w + Az) * sz + 0.5 * Ax * sx;
%% 
trf = pi / 2 / rabi;
t = 10 * pi / w;
wrf = w1;
b = [1 0 0];

v0 = U0rf(w,rabi,wrf,trf,t,b);
v1 = U1rf(Az,Ax,w,rabi,wrf,trf,t,b);
u0 = v0*v1*v1*v0;
u1 = v1*v0*v0*v1;
[a0,angle0] = axisangle(u0);
[a1,angle1] = axisangle(u1);
%%
%adjustable parameters
t = 1e-7:1e-7:1e-6;
trf = 0.1:0.1:1;%targetangle / Rabi_f / 2;
wrf = linspace(w/2,2*w1,10);

a0 = zeros(length(t),length(trf),length(wrf));
a1 = zeros(length(t),length(trf),length(wrf));
angle0 = zeros(size(axis));
angle1 = zeros(size(axis));
for i = 1:length(t)
    for j = 1:length(trf)
        for k = 1:length(wrf)
            v0 = U0rf(w,rabi,wrf(k),trf(j),t(i),b);
            v1 = U1rf(Az,Ax,w,rabi,wrf(k),trf(j),t(i),b);
            u0 = v0*v1*v1*v0;
            u1 = v1*v0*v0*v1;
            [a0(i,j,k),angle0(i,j,k)] = axisangle(u0);
            [a1(i,j,k),angle1(i,j,k)] = axisangle(u1);
        end
    end
end
ax = a0.*a1;
%% 
[a,b]=min(axis(:))
[x,y,z]=ind2sub(size(axis),b)
%% 
w = 3.3615e6;
rabi = 2 * pi * 0.1e6;
Az = 2 * pi * 10e6;
Ax = 2 * pi * 1e6;
w1 = sqrt((w + Az)^2 + Ax^2);
H0 = 0.5 * w * sz;
H1 = 0.5 * (w + Az) * sz + 0.5 * Ax * sx;

%adjustable parameters
t = 8.5e-8:0.5e-9:9.5e-8;
trf = 5e-8:0.5e-9:5.5e-8;%targetangle / Rabi_f / 2;
wrf = 6e7:0.5e6:6.5e7;

axis = zeros(length(t),length(trf),length(wrf));
angle0 = zeros(size(axis));
angle1 = zeros(size(axis));
for i = 1:length(t)
    for j = 1:length(trf)
        for k = 1:length(wrf)
            U0 = U0rf(w,rabi,wrf(k),trf(j),t(i))*U1rf(Az,Ax,w,rabi,wrf(k),2*trf(j),2*t(i))*U0rf(w,rabi,wrf(k),trf(j),t(i));
            U1 = U1rf(Az,Ax,w,rabi,wrf(k),trf(j),t(i))*U0rf(w,rabi,wrf(k),2*trf(j),2*t(i))*U1rf(Az,Ax,w,rabi,wrf(k),trf(j),t(i));
            [a0,angle0(i,j,k)]=axisangle(U0);
            [a1,angle1(i,j,k)]=axisangle(U1);
            axis(i,j,k) = a0'*a1;
        end
    end
end

%% 
%adjustable parameters
tao1 = 1e-7:1e-7:1e-5;
w_rf1 = linspace(0.8*(w + Azz),1.2*(w + Azz), 100);

axis01 = zeros(length(tao1),length(w_rf1));
axis11 = zeros(size(axis01));
angle01 = zeros(size(axis));
angle11 = zeros(size(axis));
for i = 1:length(tao1)
   for k = 1:length(w_rf1)
        U0 = expm(-1i*w_rf1*4*tao1(i)*sz/2)*U0rf1(w,rabi,w_rf1(k),tao1(i))*U1rf1(Az,Ax,w,rabi,w_rf1(k),2*tao1(i))*U0rf1(w,rabi,w_rf1(k),tao1(i));
        U1 = expm(-1i*w_rf1*4*tao1(i)*sz/2)*U1rf1(Az,Ax,w,rabi,w_rf1(k),tao1(i))*U0rf1(w,rabi,w_rf1(k),2*tao1(i))*U1rf1(Az,Ax,w,rabi,w_rf1(k),tao1(i));
        [axis01(i,k),angle01(i,k)]=axisangle(U0);
        [axis11(i,k),angle11(i,k)]=axisangle(U1);  
    end
end
%% 
%adjustable parameters
tao1 = pi/w;
tao_rf1 = 0.01:0.01:1;%targetangle / Rabi_f / 2;
w_rf1 = w/2:1e6:2*w1;

axis1 = zeros(length(tao_rf1),length(w_rf1));
angle01 = zeros(size(axis1));
angle11= zeros(size(axis1));
for j = 1:length(tao_rf1)
    for k = 1:length(w_rf1)
        U0 = U0rf(w,rabi,w_rf1(k),tao_rf1(j)*tao1,tao1)*U1rf(Az,Ax,w,rabi,w_rf1(k),2*tao_rf1(j)*tao1,2*tao1)*U0rf(w,rabi,w_rf1(k),tao_rf1(j)*tao1,tao1);
        U1 = U1rf(Az,Ax,w,rabi,w_rf1(k),tao_rf1(j)*tao1,tao1)*U0rf(w,rabi,w_rf1(k),2*tao_rf1(j)*tao1,2*tao1)*U1rf(Az,Ax,w,rabi,w_rf1(k),tao_rf1(j)*tao1,tao1);
        [a0,angle01(j,k)]=axisangle(U0);
        [a1,angle11(j,k)]=axisangle(U1);
        axis1(j,k) = a0'*a1;
    end
end
%%
[x,y]=meshgrid(tao_rf1,w_rf1);figure;surface(x,y,axis1')
[a,b]=min(axis1(:));
[x,y]=ind2sub(size(axis1),b)
%% 
%adjustable parameters
taoe = 4.45e-8:0.01e-9:4.55e-8;
w_rfe = 1.32e8:0.001e8:1.4e8;

axise = zeros(length(taoe),length(w_rfe));
angle0e = zeros(size(axise));
angle1e = zeros(size(axise));
for i = 1:length(taoe)
   for k = 1:length(w_rfe)
        U0 = U0rf1(w,rabi,w_rfe(k),taoe(i))*U1rf1(Az,Ax,w,rabi,w_rfe(k),2*taoe(i))*U0rf1(w,rabi,w_rfe(k),taoe(i));
        U1 = U1rf1(Az,Ax,w,rabi,w_rfe(k),taoe(i))*U0rf1(w,rabi,w_rfe(k),2*taoe(i))*U1rf1(Az,Ax,w,rabi,w_rfe(k),taoe(i));
        [a0,angle0e(i,k)]=axisangle(U0);
        [a1,angle1e(i,k)]=axisangle(U1);
        axise(i,k) = a0'*a1;   
    end
end
%%
[x,y]=meshgrid(taoe,w_rfe);figure;surface(x,y,axise')
[a,b]=min(axise(:));
[x,y]=ind2sub(size(axise),b);
taoe(x)
w_rfe(y)