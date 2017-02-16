w = 505*2*pi*0.00107e6;
Rabi_f = 2 * pi * 0.13e6;
A = 2 * pi * 13.7e6;
th1 = 47/180*pi;
Az = A * cos(th1);
Ax = A * sin(th1);

tao1 = linspace(0.9*(pi/A),1.1*(pi/A),100);
w_rf1 = linspace(0.9*(w+A),1.1*(w+A), 100);

axis01 = zeros(length(tao1),length(w_rf1),3);
axis11 = zeros(size(axis01));
angle01 = zeros(size(axis));
angle11 = zeros(size(axis));
for i = 1:length(tao1)
   for k = 1:length(w_rf1)
       U0 = trf01(w_rf1(k),Az,Ax,w,4*tao1(i))*U0rf1(w,Rabi_f,w_rf1(k),tao1(i))*trf10(w_rf1(k),Az,Ax,w,3*tao1(i))*U1rf1(Az,Ax,w,Rabi_f,w_rf1(k),2*tao1(i))*trf01(w_rf1(k),Az,Ax,w,tao1(i))*U0rf1(w,Rabi_f,w_rf1(k),tao1(i));
       U1 = U1rf1(Az,Ax,w,Rabi_f,w_rf1(k),tao1(i))*trf01(w_rf1(k),Az,Ax,w,3*tao1(i))*U0rf1(w,Rabi_f,w_rf1(k),2*tao1(i))*trf10(w_rf1(k),Az,Ax,w,tao1(i))*U1rf1(Az,Ax,w,Rabi_f,w_rf1(k),tao1(i));
       [axis01(i,k,:),angle01(i,k)]=axisangle(U0);
       [axis11(i,k,:),angle11(i,k)]=axisangle(U1);  
    end
end
%% 
[x,y] = meshgrid(w_rf1,tao1);
figure;surface(x,y,axis01(:,:,1)*cos(th1)-axis01(:,:,3)*sin(th1))
figure;surface(x,y,axis01(:,:,1).*axis11(:,:,1)+axis01(:,:,2).*axis11(:,:,2)+axis01(:,:,3).*axis11(:,:,3))
figure;surface(x,y,abs(angle01-angle11))
%% 
ad=axis01(:,:,1).*axis11(:,:,1)+axis01(:,:,2).*axis11(:,:,2)+axis01(:,:,3).*axis11(:,:,3);
[a,b]=min(ad(:));
[x,y]=ind2sub(size(angle01),b)
