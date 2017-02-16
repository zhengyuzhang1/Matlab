azx = sample(5,:);
t = 1e-9:1e-9:2e-5;
v0 = zeros(2,2,length(t));
v1 = zeros(size(v0));
v0dt = V0(w,t(2)-t(1));
v1dt = V1(w,t(2)-t(1),azx);
v0(:,:,1) = V0(w,t(1));
v1(:,:,1) = V1(w,t(1),azx);
for i = 2:length(t)
    v0(:,:,i) = v0dt * v0(:,:,i-1);
    v1(:,:,i) = v1dt * v1(:,:,i-1);
end
%%
t1 = 1e-8:1e-8:1e-5;
axis = zeros(length(t1),length(t1),2,3);
angle = zeros(length(t1),length(t1),2);
dotaxis = zeros(length(t1),length(t1));
idx = arrayfun(@(x)find(abs(t-x)<1e-10,1),t1);
for i = 1:length(idx)
    for j = 1:length(idx)
        tao1 = idx(i);
        tao2 = idx(j);
        u0 = v0(:,:,tao2) * v1(:,:,tao1+tao2) * v0(:,:,tao1);
        u1 = v1(:,:,tao2) * v0(:,:,tao1+tao2) * v1(:,:,tao1);
        [axis(i,j,1,:),angle(i,j,1)] = axisangle(u0);
        [axis(i,j,2,:),angle(i,j,2)] = axisangle(u1);
        dotaxis(i,j) = dot(axis(i,j,1,:),axis(i,j,2,:));
    end
end
%% 
[x,y] = meshgrid(t1(1:end/2),t1(1:end/2));
z = dotaxis(1:end/2,1:end/2);
figure;scatter(x(:),y(:),1,z(:));
%%
figure;plot(diag(dotaxis));
figure;plot(dotaxis(end/2,:));