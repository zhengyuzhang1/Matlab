n = 50;
t = linspace(1e-7,5e-6,n);
axesdot = zeros(n);
difangle = zeros(n);
condfid = zeros(n);
uncondfid = zeros(n);
for i = 1:n
    for j = 1:n
        [u0,u1,axes,angle] = generalU(w,w1{1,1}(2,:),[t(i),t(j)],1);
        axesdot(i,j) = dot(axes(:,:,1),axes(:,:,2));
        difangle(i,j) = abs(diff(angle));
        condfid(i,j) = abs(trace(u0*u1))/2;   
        uncondfid(i,j) = abs(trace(u0*u1'))/2;
    end
end
%% 
[x,y]=meshgrid(t,t);
% figure;axis square;surface(x,y,axesdot,'EdgeColor','none');
% figure;axis square;surface(x,y,difangle,'EdgeColor','none');
figure;axis square;surface(x,y,condfid,'EdgeColor','none');
figure;axis square;surface(x,y,uncondfid,'EdgeColor','none');
%% 
% figure;plot(t,diag(axesdot));
% figure;plot(t,diag(difangle));
