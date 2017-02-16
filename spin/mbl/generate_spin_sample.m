function [sortdr,indx,A,r]=generate_spin_sample(N,smallestdistance)
r=zeros(N,3);
dx=zeros(N,N,3);
dr=zeros(N,N);
A=zeros(N,N,6);
n=0;
while n<N
   n=n+1;
   r(n,:)=rand(1,3);
   for i=1:n-1
       dx(i,n,:)=r(i,:)-r(n,:);
       dr(i,n)=sqrt(dx(i,n,1)^2+dx(i,n,2)^2+dx(i,n,3)^2);
       if dr(i,n)<smallestdistance
           n=n-1;
           break
       end
   end
end
for i=1:N
    for j=i+1:N
        dx(i,j,:)=dx(i,j,:)/dr(i,j);
        A(i,j,:)=-1/dr(i,j)*[3*dx(i,j,1)^2-1,3*dx(i,j,1)*dx(i,j,2),3*dx(i,j,1)*dx(i,j,3),3*dx(i,j,2)^2-1,3*dx(i,j,2)*dx(i,j,3),3*dx(i,j,3)^2-1];
        A(j,i,:)=A(i,j,:);    
    end
end        
[sortdr,indx]=sort(dr(:));
indx=indx(sortdr>0);
sortdr=sortdr(sortdr>0);
mask=sortdr<median(sortdr);
sortdr=sortdr(mask);
indx=indx(mask,:);
[indI,indJ]=ind2sub([N,N],indx);
indx=[indI,indJ];
end